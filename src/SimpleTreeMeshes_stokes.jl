using SparseArrays
using LinearAlgebra

@enum StokesBoundaryConditions stokes_boundary_free_slip stokes_boundary_dirichlet

function sol2vp(tree::SimpleTreeMesh, x::Vector{Float64}, Kcont::Float64)
    @assert length(x) == tree.nf + tree.ne
    return x[1:tree.nf], x[tree.nf+1:end] * Kcont
end;

function vp2sol(tree::SimpleTreeMesh, v::Vector{Float64}, p::Vector{Float64}, Kcont::Float64)
    @assert length(v) == tree.nf
    @assert length(p) == tree.ne
    return [v; p/Kcont]
end

function e_id_to_eqnum(tree::SimpleTreeMesh, e_id::Int)
     @assert 1 <= e_id <= tree.ne
     return e_id + tree.nf
end;

function e_id_to_eqnum(tree::SimpleTreeMesh, e_id_array::Array{Int})
    return map(e_id -> e_id_to_eqnum(tree, e_id), e_id_array)
end;

function f_id_to_eqnum(tree::SimpleTreeMesh, f_id::Int)
     @assert 1 <= f_id <= tree.nf
     return f_id
end;

function f_id_to_eqnum(tree::SimpleTreeMesh, f_id_array::Array{Int})
    return map(f_id -> f_id_to_eqnum(tree, f_id), f_id_array)
end;

function eqnum_to_f_id(tree::SimpleTreeMesh, eqnum::Int)
     @assert 1 <= eqnum <= tree.nf
     return eqnum
end;

function eqnum_to_e_id(tree::SimpleTreeMesh, eqnum::Int)
     @assert tree.nf + 1 <= eqnum <= tree.nf + tree.ne
     return eqnum - tree.nf
end;


function dvx_dy(tree::SimpleTreeMesh, c_id, problem=nothing)
    @assert 1 <= c_id <= tree.nc
    corner_type = tree.c[c_id].corner_type
    corner_info = tree.c[c_id].info

    if corner_type == corner_boundary
        @assert problem.boundary_conditions in [stokes_boundary_free_slip, stokes_boundary_dirichlet]
        if corner_info & 0b1100 == 0b1100  # One of the corners of the domain
            # This domain corner case is not required for free-slip or Dirichlet BCs,
            # just when wanting to directly compute Txy (e.g. for testing)
            vx_boundary =  problem.vx_ref(get_corner_coordinates(tree, c_id)...)
            if corner_info in (CORNER_BOUNDARY_DOWN_LEFT, CORNER_BOUNDARY_DOWN_RIGHT)
                f_id_up = tree.c2f[UP, c_id]
                h_up = get_face_size(tree, f_id_up) / 2
                cols = f_id_to_eqnum(tree, [f_id_up])
                vals =  [1.0/h_up]
                rhs_val = vx_boundary / h_up
            else
                @assert corner_info in (CORNER_BOUNDARY_UP_RIGHT, CORNER_BOUNDARY_UP_LEFT)
                f_id_down = tree.c2f[DOWN, c_id]
                h_down = get_face_size(tree, f_id_down) / 2
                cols = f_id_to_eqnum(tree, [f_id_down])
                vals = [-1.0/h_down]
                rhs_val = -1.0 * vx_boundary / h_down
            end
        else
            # On the boundaries, we always make the choice to use small faces, when split
            corner_type_masked = corner_info & 0b0011
            if corner_type_masked in [CORNER_BOUNDARY_DOWN_MASK, CORNER_BOUNDARY_UP_MASK]
                if problem.boundary_conditions == stokes_boundary_free_slip
                    # Free-slip boundary conditions have dvx_dy = 0 on the top and bottom
                    cols = []
                    vals = []
                    rhs_val = 0.0
                else
                    @assert problem.boundary_conditions == stokes_boundary_dirichlet
                    vx_boundary =  problem.vx_ref(get_corner_coordinates(tree, c_id)...)
                    if corner_type_masked == CORNER_BOUNDARY_DOWN_MASK
                        f_id_up =c2f_small(tree, UP, c_id)
                        cols = f_id_to_eqnum(tree, [f_id_up])
                        h_up = get_face_size(tree, f_id_up)/2
                        vals =  [1.0/h_up]
                        rhs_val = vx_boundary / h_up
                    else
                        @assert corner_type_masked == CORNER_BOUNDARY_UP_MASK
                        f_id_down = c2f_small(tree, DOWN, c_id)
                        cols = f_id_to_eqnum(tree, [f_id_down])
                        h_down = get_face_size(tree, f_id_down)/2
                        vals = [-1.0/h_down]
                        rhs_val = -1.0 * vx_boundary / h_down
                    end
                end
            else
                @assert corner_type_masked in (CORNER_BOUNDARY_LEFT_MASK, CORNER_BOUNDARY_RIGHT_MASK)
                f_id_up = c2f_small(tree, UP, c_id)  # should actually never be small
                h_up = get_face_size(tree, f_id_up)/2
                f_id_down = c2f_small(tree, DOWN, c_id)  # should actually never be small
                h_down = get_face_size(tree, f_id_down)/2
                cols = f_id_to_eqnum(tree, [f_id_down, f_id_up])
                hy = h_down + h_up
                vals = [-1.0/hy, 1.0/hy]
                rhs_val = 0.0
            end
        end
    elseif corner_type == corner_regular
        # following GMD13, in most cases we just take the closest velocities,
        # except in the case where we have exactly 3 coarse cells and one fine cell
        choose_coarse = corner_info in [0b0001, 0b0010, 0b0100, 0b1000]

        # for convenience, split info on the level (0=coarse, 1=fine) of the 4 neighboring elements
        down_right_is_fine = Bool((corner_info & 0b0100) >> 2)
        up_right_is_fine   = Bool((corner_info & 0b0001) >> 0)

        if xor(choose_coarse, down_right_is_fine)
            e_id_down = tree.c2e[DOWN_RIGHT, c_id]
            f_id_down = tree.e2f[LEFT, e_id_down]
        else
            e_id_down = tree.c2e[DOWN_LEFT, c_id]
            f_id_down = tree.e2f[RIGHT, e_id_down]
        end
        h_down = get_element_size(tree, e_id_down) / 2

        if xor(choose_coarse, up_right_is_fine)
            e_id_up = tree.c2e[UP_RIGHT, c_id]
            f_id_up = tree.e2f[LEFT, e_id_up]
        else
            e_id_up = tree.c2e[UP_LEFT, c_id]
            f_id_up = tree.e2f[RIGHT, e_id_up]
        end
        h_up = get_element_size(tree, e_id_up) / 2

        cols = f_id_to_eqnum(tree, [f_id_down, f_id_up])
        hy = h_down + h_up
        vals = [-1.0/hy, 1.0/hy]
        rhs_val = 0.0
    else
        @assert corner_type == corner_hanging
        # In GMD13, Edot_xy is never computed directly on hanging corners,
        # but we include this "method 2" for comparison's sake

        # "method 2" approach - on big element, get missing velocity by linear interpolation=averaging,
        # and then apply a local viscosity. Interestingly this seems to be less accurate than "method 1"
        # as in GMD13

        e_id_big = hanging_corner_get_big_element(tree, c_id)

        if corner_info == CORNER_HANG_LEFT
            e_id_up_right = tree.c2e[UP_RIGHT, c_id]
            e_id_down_right = tree.c2e[DOWN_RIGHT, c_id]
            h_fine = get_element_size(tree, e_id_up_right) / 2
            h_up = h_down = h_fine
            f_ids = [
                tree.e2f[LEFT, e_id_down_right],
                tree.e2f[LEFT, e_id_up_right],
            ]
            hy = h_down + h_up
            vals = [ -1.0/hy, 1.0/hy]
            rhs_val = 0.0
        elseif corner_info == CORNER_HANG_RIGHT
            e_id_up_left = tree.c2e[UP_LEFT, c_id]
            e_id_down_left = tree.c2e[DOWN_LEFT, c_id]
            h_fine = get_element_size(tree, e_id_up_left) / 2
            h_up = h_down = h_fine
            f_ids = [
                tree.e2f[RIGHT, e_id_down_left],
                tree.e2f[RIGHT, e_id_up_left],
            ]
            hy = h_down + h_up
            vals = [-1.0/hy, 1.0/hy]
            rhs_val = 0.0
        elseif corner_info == CORNER_HANG_DOWN
            e_id_up_right = tree.c2e[UP_RIGHT, c_id]
            e_id_up_left = tree.c2e[UP_LEFT, c_id]
            h_fine = get_element_size(tree, e_id_up_right) / 2
            h_up = h_fine
            h_down = 2 * h_fine
            f_ids = [
                tree.e2f[LEFT, e_id_big],
                tree.e2f[RIGHT, e_id_big],
                tree.e2f[LEFT, e_id_up_right],
            ]
            hy = h_down + h_up
            vals = [-1.0/(2*hy), -1.0/(2*hy), 1.0/hy]
            rhs_val = 0.0
        else
            @assert corner_info == CORNER_HANG_UP
            e_id_down_left = tree.c2e[DOWN_LEFT, c_id]
            e_id_down_right = tree.c2e[DOWN_RIGHT, c_id]
            h_fine = get_element_size(tree, e_id_down_left) / 2
            h_down = h_fine
            h_up = 2 * h_fine
            f_ids = [
                tree.e2f[LEFT, e_id_down_right],
                tree.e2f[LEFT, e_id_big],
                tree.e2f[RIGHT, e_id_big],
            ]
            hy = h_down + h_up
            vals = [-1.0/hy, 1.0/(2*hy), 1.0/(2*hy)]
            rhs_val = 0.0
        end
        cols = f_id_to_eqnum(tree, f_ids)
    end
    return cols, vals, rhs_val
end;

function dvy_dx(tree::SimpleTreeMesh, c_id, problem=nothing)
    @assert 1 <= c_id <= tree.nc
    corner_type = tree.c[c_id].corner_type
    corner_info = tree.c[c_id].info

    if corner_type == corner_boundary
        @assert problem.boundary_conditions in [stokes_boundary_free_slip, stokes_boundary_dirichlet]
        if corner_info & 0b1100 == 0b1100
            # This domain corner case is not actually required for free-slip or Dirichlet BCs
            vy_boundary =  problem.vy_ref(get_corner_coordinates(tree, c_id)...)
            if corner_info in (CORNER_BOUNDARY_DOWN_LEFT, CORNER_BOUNDARY_UP_LEFT)
                f_id_right = tree.c2f[RIGHT, c_id]
                h_right = get_face_size(tree, f_id_right) / 2
                cols = f_id_to_eqnum(tree, [f_id_right])
                vals = [1.0/h_right]
                rhs_val = vy_boundary / h_right
            else
                @assert corner_info in (CORNER_BOUNDARY_UP_RIGHT, CORNER_BOUNDARY_DOWN_RIGHT)
                f_id_left = tree.c2f[LEFT, c_id]
                h_left = get_face_size(tree, f_id_left) / 2
                cols = f_id_to_eqnum(tree, [f_id_left])
                vals = -1.0 / h_left
                rhs_val = -1.0 * vy_boundary / h_left
            end
        else
            corner_type_masked = corner_info & 0b0011
            if corner_type_masked in [CORNER_BOUNDARY_LEFT_MASK, CORNER_BOUNDARY_RIGHT_MASK]
                if problem.boundary_conditions == stokes_boundary_free_slip
                    # Free-slip boundary conditions with dvy_dx = 0 on the left and right
                    cols = []
                    vals = []
                    rhs_val = 0.0
                else
                    @assert problem.boundary_conditions == stokes_boundary_dirichlet
                    vy_boundary =  problem.vy_ref(get_corner_coordinates(tree, c_id)...)
                    if corner_type_masked == CORNER_BOUNDARY_LEFT_MASK
                        f_id_right = c2f_small(tree, RIGHT, c_id)
                        cols = f_id_to_eqnum(tree, [f_id_right])
                        h_right = get_face_size(tree, f_id_right)/2
                        vals = [1.0/h_right]
                        rhs_val = vy_boundary / h_right
                    else
                        @assert corner_type_masked == CORNER_BOUNDARY_RIGHT_MASK
                        f_id_left = c2f_small(tree, LEFT, c_id)
                        cols = f_id_to_eqnum(tree, [f_id_left])
                        h_left = get_face_size(tree, f_id_left)/2
                        vals = [-1.0/h_left]
                        rhs_val = -1.0 * vy_boundary / h_left
                    end
                end
            else
                @assert corner_type_masked in (CORNER_BOUNDARY_DOWN_MASK, CORNER_BOUNDARY_UP_MASK)
                f_id_left = c2f_small(tree, LEFT, c_id)  # should actually never be small
                h_left = get_face_size(tree, f_id_left) / 2
                f_id_right = c2f_small(tree, RIGHT, c_id)  # should actually never be small
                h_right = get_face_size(tree, f_id_right) / 2
                cols = f_id_to_eqnum(tree, [f_id_left, f_id_right])
                hx = h_left + h_right
                vals = [-1.0/hx, 1.0/hx]
                rhs_val = 0.0
            end
        end
    elseif corner_type == corner_regular
        choose_coarse = corner_info in [0b0001, 0b0010, 0b0100, 0b1000]

        up_left_is_fine    = Bool((corner_info & 0b0010) >> 1)
        up_right_is_fine   = Bool((corner_info & 0b0001) >> 0)

        if xor(choose_coarse, up_left_is_fine) # if exactly one is true
            e_id_left = tree.c2e[UP_LEFT, c_id]
            f_id_left = tree.e2f[DOWN, e_id_left]
        else
            e_id_left= tree.c2e[DOWN_LEFT, c_id]
            f_id_left = tree.e2f[UP, e_id_left]
        end
        h_left = get_element_size(tree, e_id_left) / 2

        if xor(choose_coarse, up_right_is_fine)
            e_id_right = tree.c2e[UP_RIGHT, c_id]
            f_id_right = tree.e2f[DOWN, e_id_right]
        else
            e_id_right = tree.c2e[DOWN_RIGHT, c_id]
            f_id_right = tree.e2f[UP, e_id_right]
        end
        h_right = get_element_size(tree, e_id_right) / 2

        cols = f_id_to_eqnum(tree, [f_id_left, f_id_right])
        hx = h_left + h_right
        vals = [-1.0/hx, 1.0/hx]
        rhs_val = 0.0
    else
        @assert corner_type == corner_hanging
        # "method 2" approach - on big element, get missing velocity by linear interpolation=averaging,
        # and then apply a local viscosity. Interestingly this seems to be less accurate than "method 1"
        # as in GMD13

        e_id_big = hanging_corner_get_big_element(tree, c_id)

        if corner_info == CORNER_HANG_LEFT
            e_id_up_right = tree.c2e[UP_RIGHT, c_id]
            e_id_down_right = tree.c2e[DOWN_RIGHT, c_id]
            h_fine = get_element_size(tree, e_id_up_right) / 2
            h_right = h_fine
            h_left = 2 * h_fine
            f_ids = [
                tree.e2f[UP, e_id_big],
                tree.e2f[DOWN, e_id_big],
                tree.e2f[DOWN, e_id_up_right],
            ]
            hx = h_left + h_right
            vals = [-1.0/(2*hx), -1.0/(2*hx), 1.0/hx]
            rhs_val = 0.0
        elseif corner_info == CORNER_HANG_RIGHT
            e_id_up_left = tree.c2e[UP_LEFT, c_id]
            e_id_down_left = tree.c2e[DOWN_LEFT, c_id]
            h_fine = get_element_size(tree, e_id_up_left) / 2
            h_left = h_fine
            h_right = 2 * h_fine
            f_ids = [
                tree.e2f[DOWN, e_id_up_left],
                tree.e2f[UP, e_id_big],
                tree.e2f[DOWN, e_id_big],
            ]
            hx = h_left + h_right
            vals = [-1.0/hx, 1.0/(2*hx), 1.0/(2*hx)]
            rhs_val = 0.0
        elseif corner_info == CORNER_HANG_DOWN
            e_id_up_right = tree.c2e[UP_RIGHT, c_id]
            e_id_up_left = tree.c2e[UP_LEFT, c_id]
            h_fine = get_element_size(tree, e_id_up_right) / 2
            h_left = h_right = h_fine
            f_ids = [
                tree.e2f[DOWN, e_id_up_left],
                tree.e2f[DOWN, e_id_up_right],
            ]
            hx = h_left + h_right
            vals = [-1.0/hx, 1.0/hx]
            rhs_val = 0.0
        else
            @assert corner_info == CORNER_HANG_UP
            e_id_down_left = tree.c2e[DOWN_LEFT, c_id]
            e_id_down_right = tree.c2e[DOWN_RIGHT, c_id]
            h_fine = get_element_size(tree, e_id_down_left) / 2
            h_left = h_right = h_fine
            f_ids = [
                tree.e2f[UP, e_id_down_left],
                tree.e2f[UP, e_id_down_right],
            ]
            hx = h_left + h_right
            vals = [-1.0/hx, 1.0/hx]
            rhs_val = 0.0
        end
        cols = f_id_to_eqnum(tree, f_ids)
    end
    return cols, vals, rhs_val
end;

function two_times_Edot_xy_terms(tree::SimpleTreeMesh, c_id, problem=nothing)
    @assert 1 <= c_id <= tree.nc
    cols1, vals1, rhs_val1 = dvx_dy(tree, c_id, problem)
    cols2, vals2, rhs_val2 = dvy_dx(tree, c_id, problem)
    return [cols1; cols2], [vals1; vals2], rhs_val1 + rhs_val2
end;

function compute_Txy(tree::SimpleTreeMesh, problem, v::Vector{Float64})
  @assert length(v) == tree.nf

  Txy = fill(Float64(0), tree.nc)
  for c_id = 1:tree.nc
    corner_type = tree.c[c_id].corner_type

    if corner_type == corner_hanging
      error("Not implemented for hanging corners")
      # Note: if implementing this, note different "Txy_method" values used in assemble_system
    else
      @assert corner_type == corner_boundary || corner_type == corner_regular
      cols, vals, rhs_val = two_times_Edot_xy_terms(tree, c_id, problem)
      scale = problem.eta(get_corner_coordinates(tree, c_id)...)
      for i = 1:length(cols)
        f_id = cols[i]  # cols is for a v-p vector
        Txy[c_id] += scale * vals[i] * v[f_id]
      end
      Txy[c_id] -= scale * rhs_val
    end
  end
  return Txy
end

# Note: the underlying p4est structure can tell you h_min, so that argument could be eliminated with some work
function assemble_system(tree::SimpleTreeMesh, problem, h_min; Txy_method=1, zero_diri_vel_cols=false, nonzero_x0=false, modified_scaling=false, modified_stress=false)

    # Check for a single boundary condition type in the problem
    # (future development could specify this on a per-boundary or even per-cell basis)
    @assert problem.boundary_conditions in [stokes_boundary_dirichlet, stokes_boundary_free_slip]

    n = tree.nf + tree.ne
    row = Vector{Int64}()
    col = Vector{Int64}()
    val = Vector{Float64}()
    b = Vector{Float64}(undef, n)

    if modified_scaling
      Kcont = 1.0
      Kbnd = 1.0
    else
      Kcont = problem.eta_ref / h_min
      Kbnd =  problem.eta_ref / (h_min^2)
    end

    for e_id = 1:tree.ne
        eqnum = e_id_to_eqnum(tree, e_id)
        if e_id == 1
            # Fix/pin single pressure node
            push!(row, eqnum)
            push!(col, eqnum)
            push!(val, Kbnd)
            b[eqnum] = (Kbnd / Kcont) * problem.p_ref(get_element_coordinates(tree, e_id)...)
        else
            # Mass conservation terms
            push!(row, eqnum, eqnum, eqnum, eqnum)
            push!(col, f_id_to_eqnum(tree, tree.e2f[:, e_id])...)
            hx = hy = get_element_size(tree, e_id)
            push!(val, -Kcont/hx, Kcont/hx, -Kcont/hy, Kcont/hy)
            b[eqnum] = Kcont * problem.f_p(get_element_coordinates(tree, e_id)...)
        end

    end

    for f_id = 1:tree.nf
        eqnum = f_id_to_eqnum(tree, f_id)
        b[eqnum] = 0.0
        face_type = tree.f[f_id].face_type
        face_direction = tree.f[f_id].direction
        if face_type == face_boundary
            # Dirichlet boundaries
            push!(row, eqnum)
            push!(col, eqnum)
            push!(val, Kbnd)
            if face_direction == face_vx
                b[eqnum] = Kbnd * problem.vx_ref(get_face_coordinates(tree, f_id)...)
            else
                @assert face_direction == face_vy
                b[eqnum] = Kbnd * problem.vy_ref(get_face_coordinates(tree, f_id)...)
            end
        elseif face_type == face_regular
            e_id_lo, e_id_hi = tree.f2e[:, f_id]

            # - Grad p terms
            push!(row, eqnum, eqnum)
            @assert !(-1 in tree.f2e[:, f_id])
            push!(col, e_id_to_eqnum(tree, [e_id_lo, e_id_hi])...)
            h = get_element_size(tree, e_id_lo)
            push!(val,  Kcont/h, -Kcont/h)

            # Div(stress) terms

            # Txx/Tyy
            h2 = h^2
            LO, HI = face_direction == face_vx ? (LEFT, RIGHT) : (DOWN, UP)
            f_id_lo = tree.e2f[LO, e_id_lo]
            f_id_hi = tree.e2f[HI, e_id_hi]
            if modified_stress
              LO2, HI2 = face_direction == face_vx ? (DOWN, UP) : (LEFT, RIGHT)
            end
            eta_lo = problem.eta(get_element_coordinates(tree, e_id_lo)...)
            eta_hi = problem.eta(get_element_coordinates(tree, e_id_hi)...)


            if modified_stress
              f_id_lo2 = tree.e2f[LO2, e_id_lo] # both for lo
              f_id_hi2 = tree.e2f[HI2, e_id_lo] # both for lo
              push!(row, eqnum, eqnum, eqnum, eqnum)
              push!(col, f_id_to_eqnum(tree, [f_id_lo, f_id, f_id_lo2, f_id_hi2])...)
              push!(val, (4.0/3.0)*eta_lo/h2, -(4.0/3.0)*eta_lo/h2, -(2.0/3.0)*eta_lo/h2, (2.0/3.0)*eta_lo/h2)
              # Note: h2 is used for all of hx^2, hy^2, hxhy, assuming square elements
            else
              push!(row, eqnum, eqnum)
              push!(col, f_id_to_eqnum(tree, [f_id_lo, f_id])...)
              push!(val, 2.0*eta_lo/h2, -2.0*eta_lo/h2)
            end

            if modified_stress
              f_id_lo2 = tree.e2f[LO2, e_id_hi] # both for hi
              f_id_hi2 = tree.e2f[HI2, e_id_hi] # both for hi
              push!(row, eqnum, eqnum, eqnum, eqnum)
              push!(col, f_id_to_eqnum(tree, [f_id, f_id_hi, f_id_lo2, f_id_hi2])...)
              push!(val, -(4.0/3.0)*eta_hi/h2, (4.0/3.0)*eta_hi/h2, (2.0/3.0)*eta_hi/h2, -(2.0/3.0)*eta_hi/h2)
              # Note: h2 is used for all of hx^2, hy^2, hxhy, assuming square elements
            else
              push!(row, eqnum, eqnum)
              push!(col, f_id_to_eqnum(tree, [f_id, f_id_hi])...)
              push!(val, -2.0*eta_hi/h2, 2.0*eta_hi/h2)
            end

            # Txy

            # FIXME: redundant logic for lo2 and hi2
            c_id_lo2 = tree.e2c[DOWN_LEFT, e_id_hi]
            if Txy_method !=2 && tree.c[c_id_lo2].corner_type == corner_hanging
                @assert tree.c[c_id_lo2].info == (face_direction == face_vx ? CORNER_HANG_DOWN : CORNER_HANG_LEFT)
                # Compute Txy on hanging corner by averaging Txy on adjacent non-hanging corners
                e_id_big = hanging_corner_get_big_element(tree, c_id_lo2)
                c_positions = face_direction == face_vx ?  (UP_LEFT, UP_RIGHT) : (UP_RIGHT, DOWN_RIGHT)
                if Txy_method == 3
                  eta_lo2 = problem.eta(get_corner_coordinates(tree, c_id_lo2)...)
                end
                for c_position in c_positions
                    c_id_contrib = tree.e2c[c_position, e_id_big]
                    cols_contrib, vals_contrib, rhs_val_contrib = two_times_Edot_xy_terms(tree, c_id_contrib, problem)
                    if Txy_method == 3
                      # "method 3" means apply a local viscosity
                      eta_contrib = eta_lo2
                    else
                      eta_contrib = problem.eta(get_corner_coordinates(tree, c_id_contrib)...)
                    end
                    push!(row, fill(eqnum, length(cols_contrib))...)
                    push!(col, cols_contrib...)
                    val_scale = -eta_contrib / (2*h)
                    push!(val, (vals_contrib * val_scale)...)
                    b[eqnum] += rhs_val_contrib * val_scale
                end
            else
                cols_lo2, vals_lo2, rhs_val_lo2 = two_times_Edot_xy_terms(tree, c_id_lo2, problem)
                eta_lo2 = problem.eta(get_corner_coordinates(tree, c_id_lo2)...)
                push!(row, fill(eqnum, length(cols_lo2))...)
                push!(col, cols_lo2...)
                val_scale = -eta_lo2 / h
                push!(val, (vals_lo2 * val_scale)...)
                b[eqnum] += rhs_val_lo2 * val_scale
            end

            c_id_hi2 = tree.e2c[UP_RIGHT, e_id_lo]
            if Txy_method != 2 && tree.c[c_id_hi2].corner_type == corner_hanging
                @assert tree.c[c_id_hi2].info == (face_direction == face_vx ? CORNER_HANG_UP : CORNER_HANG_RIGHT)
                e_id_big = hanging_corner_get_big_element(tree, c_id_hi2)
                c_positions = face_direction == face_vx ?  (DOWN_LEFT, DOWN_RIGHT) : (UP_LEFT, DOWN_LEFT)
                if Txy_method == 3
                  eta_hi2 = problem.eta(get_corner_coordinates(tree, c_id_hi2)...)
                end
                for c_position in c_positions
                    c_id_contrib = tree.e2c[c_position, e_id_big]
                    cols_contrib, vals_contrib, rhs_val_contrib = two_times_Edot_xy_terms(tree, c_id_contrib, problem)
                    if Txy_method == 3
                      eta_contrib = eta_hi2
                    else
                      eta_contrib = problem.eta(get_corner_coordinates(tree, c_id_contrib)...)
                    end
                    push!(row, fill(eqnum, length(cols_contrib))...)
                    push!(col, cols_contrib...)
                    val_scale = eta_contrib/(2*h)
                    push!(val, (vals_contrib * val_scale)...)
                    b[eqnum] += rhs_val_contrib * val_scale
                end
            else
                cols_hi2, vals_hi2, rhs_val_hi2 = two_times_Edot_xy_terms(tree, c_id_hi2, problem)
                eta_hi2 = problem.eta(get_corner_coordinates(tree, c_id_hi2)...)
                push!(row, fill(eqnum,length(cols_hi2))...)
                push!(col, cols_hi2...)
                val_scale = eta_hi2/h
                push!(val, (vals_hi2 * val_scale)...)
                b[eqnum] += rhs_val_hi2 * val_scale
            end

            if face_direction == face_vx
                b[eqnum] += problem.f_x(get_face_coordinates(tree, f_id)...)
            else
                @assert face_direction == face_vy
                b[eqnum] += problem.f_y(get_face_coordinates(tree, f_id)...)
            end
        elseif face_type == face_big
            big_low = tree.f2e[2, f_id] == -1 # The non-(-1) entry is first
            big_slot = big_low ? 1 : 2
            small_slot = big_low ? 2 : 1
            e_id_big = tree.f2e[big_slot, f_id]
            e_id_first = tree.f2e[small_slot, f_id+1] # brittle
            e_id_second = tree.f2e[small_slot, f_id+2] # brittle

            # - Grad P terms
            push!(row, eqnum, eqnum, eqnum)
            push!(col, e_id_to_eqnum(tree, [e_id_big, e_id_first, e_id_second])...)
            h = get_element_size(tree, e_id_big) * (3/4)
            if big_low
               push!(val, Kcont/h, -0.5*Kcont/h, -0.5*Kcont/h)
            else
               push!(val, -Kcont/h, 0.5*Kcont/h, 0.5*Kcont/h)
            end

            # Txx/Tyy terms
            if modified_stress; @error("Not implemented!"); end
            LO, HI = face_direction == face_vx ? (LEFT, RIGHT) : (DOWN, UP)
            f_id_big_lo = tree.e2f[LO, e_id_big]
            f_id_big_hi = tree.e2f[HI, e_id_big]
            @assert (big_low && f_id_big_hi == f_id) || f_id_big_lo == f_id
            f_id_first_lo = tree.e2f[LO, e_id_first]
            f_id_first_hi = tree.e2f[HI, e_id_first]
            f_id_second_lo = tree.e2f[LO, e_id_second]
            f_id_second_hi = tree.e2f[HI, e_id_second]

            push!(row, fill(eqnum, 6)...)
            push!(col, f_id_to_eqnum(tree, [f_id_big_lo, f_id_big_hi, f_id_first_lo, f_id_first_hi, f_id_second_lo, f_id_second_hi])...)
            h = get_element_size(tree, e_id_big) * (3/4) # redundant
            h_big = get_element_size(tree, e_id_big)
            h_small = get_element_size(tree, e_id_first)
            @assert get_element_size(tree, e_id_second) == h_small
            eta_big = problem.eta(get_element_coordinates(tree, e_id_big)...)
            eta_first = problem.eta(get_element_coordinates(tree, e_id_first)...)
            eta_second = problem.eta(get_element_coordinates(tree, e_id_second)...)
            s = big_low ? 1.0 : -1.0
            push!(val, s * 2.0 * eta_big / (h*h_big), -s*2.0*eta_big/(h*h_big), -s*eta_first / (h*h_small), s*eta_first/(h*h_small), -s*eta_second/(h*h_small), s*eta_second/(h*h_small))

            # Txy terms
            if face_direction == face_vx
                if big_low
                    c_id_lo2, c_id_hi2 = tree.e2c[[DOWN_RIGHT, UP_RIGHT], e_id_big]
                else
                    c_id_lo2, c_id_hi2 = tree.e2c[[DOWN_LEFT, UP_LEFT], e_id_big]
                end
            else
                @assert face_direction == face_vy
                if big_low
                    c_id_lo2, c_id_hi2 = tree.e2c[[UP_LEFT, UP_RIGHT], e_id_big]
                else
                    c_id_lo2, c_id_hi2 = tree.e2c[[DOWN_LEFT, DOWN_RIGHT], e_id_big]
                end
            end
            @assert tree.c[c_id_lo2].corner_type != corner_hanging
            @assert tree.c[c_id_hi2].corner_type != corner_hanging

            cols_lo2, vals_lo2, rhs_val_lo2 = two_times_Edot_xy_terms(tree, c_id_lo2, problem)
            eta_lo2 = problem.eta(get_corner_coordinates(tree, c_id_lo2)...)
            h = get_element_size(tree, e_id_big)

            push!(row, fill(eqnum, length(cols_lo2))...)
            push!(col, cols_lo2...)
            val_scale = -eta_lo2/h
            push!(val, (vals_lo2 * val_scale)...)
            b[eqnum] += rhs_val_lo2 * val_scale

            cols_hi2, vals_hi2, rhs_val_hi2 = two_times_Edot_xy_terms(tree, c_id_hi2, problem)
            eta_hi2 = problem.eta(get_corner_coordinates(tree, c_id_hi2)...)
            push!(row, fill(eqnum,length(cols_hi2))...)
            push!(col, cols_hi2...)
            val_scale = eta_hi2/h
            push!(val, (vals_hi2 * val_scale)...)
            b[eqnum] += rhs_val_hi2 * val_scale

            if face_direction == face_vx
                b[eqnum] += problem.f_x(get_face_coordinates(tree, f_id)...)
            else
                @assert face_direction == face_vy
                b[eqnum] += problem.f_y(get_face_coordinates(tree, f_id)...)
            end
        else
            @assert face_type == face_first || face_type == face_second
            if modified_stress; @error("Not implemented!"); end  # haven't though about it, anyway
            # FIXME there is probaby still too much redundancy with the "big" case
            if face_type == face_first
                f_id_big = f_id-1
            else
                @assert face_type == face_second
                f_id_big = f_id-2
            end
            big_low = tree.f2e[2, f_id_big] == -1
            e_id_big = big_low ? tree.f2e[1, f_id_big] : tree.f2e[2, f_id_big]

            if face_direction == face_vx
                if big_low
                    c_id_lo2, c_id_hi2 = tree.e2c[[DOWN_RIGHT, UP_RIGHT], e_id_big]
                else
                    c_id_lo2, c_id_hi2 = tree.e2c[[DOWN_LEFT, UP_LEFT], e_id_big]
                end
            else
                @assert face_direction == face_vy
                if big_low
                    c_id_lo2, c_id_hi2 = tree.e2c[[UP_LEFT, UP_RIGHT], e_id_big]
                else
                    c_id_lo2, c_id_hi2 = tree.e2c[[DOWN_LEFT, DOWN_RIGHT], e_id_big]
                end
            end
            @assert tree.c[c_id_lo2].corner_type != corner_hanging
            @assert tree.c[c_id_hi2].corner_type != corner_hanging

            # FIXME: We'd like to relax this, but have to deal with eliminated Diri conditions, I think!
            @assert tree.c[c_id_lo2].corner_type != corner_boundary || problem.boundary_conditions != stokes_boundary_dirichlet
            @assert tree.c[c_id_hi2].corner_type != corner_boundary || problem.boundary_conditions != stokes_boundary_dirichlet

            eta_lo2 = problem.eta(get_corner_coordinates(tree, c_id_lo2)...)
            eta_hi2 = problem.eta(get_corner_coordinates(tree, c_id_hi2)...)
            h_big = get_element_size(tree, e_id_big)

            if face_direction == face_vx
                cols_hi2, vals_hi2 = dvx_dy(tree, c_id_hi2, problem)
            else
                @assert face_direction == face_vy
                cols_hi2, vals_hi2 = dvy_dx(tree, c_id_hi2, problem)
            end
            if length(cols_hi2) > 0
                @assert length(cols_hi2) == 2
                @assert length(vals_hi2) == 2
                push!(row, eqnum, eqnum)
                push!(col, cols_hi2...)
                push!(val, (-1.0 * eta_hi2 * vals_hi2)...)
            end

            if face_direction == face_vx
                cols_lo2, vals_lo2 = dvx_dy(tree, c_id_lo2, problem)
            else
                @assert face_direction == face_vy
                cols_lo2, vals_lo2 = dvy_dx(tree, c_id_lo2, problem)
            end
            if length(cols_lo2) > 0
                @assert length(cols_lo2) == 2
                @assert length(vals_lo2) == 2
                push!(row, eqnum, eqnum)
                push!(col, cols_lo2...)
                push!(val, (-1.0 * eta_lo2 * vals_lo2)...)
            end

            eta_hat = eta_lo2 + eta_hi2  # includes factor of 2
            push!(row, eqnum, eqnum)
            push!(col, f_id_to_eqnum(tree, [f_id, f_id_big])...)
            sign = face_type == face_first ? 1.0 : -1.0
            push!(val, [-sign*eta_hat/(h_big/4), sign*eta_hat/(h_big/4)]...)

        end
    end

    A = sparse(row, col, val)

    # Form a system which solves for the update from a state with Dirichlet
    # velocity conditions already specified
    if zero_diri_vel_cols || nonzero_x0
      @assert length(col) == length(row)
      @assert length(val) == length(col)
      for i = 1:length(col)
        if col[i] <= tree.nf
          f_id = eqnum_to_f_id(tree, col[i])
          if tree.f[f_id].face_type == face_boundary
            val[i] = col[i] == row[i] ? Kbnd : 0.0
          end
        end
      end

      if nonzero_x0
        x0 = populate_x0(tree, problem)
      else
        x = fill(Float64(0), tree.nf + tree.ne)
      end
      populate_stokes_boundary_velocities!(tree, problem, x0)
      bmod = b - (A * x0)
      Amod = sparse(row, col, val)
      return Amod, bmod, Kcont, x0
    else
      return A, b, Kcont, nothing
    end

end;

function populate_x0(tree, problem)
  # Populate as in TD's SetInitialVelocity2()
  # This version is wrong on the boundary! (overwrite with boundary values later)
  # We sum two linear functions, as TD does, but maybe it'd be better to use a bilinear function
  function vy0(x::Float64, y::Float64)
    x_min = tree.coordinate_offset_x
    y_min = tree.coordinate_offset_y
    x_max = tree.coordinate_offset_x + tree.coordinate_scale
    y_max = tree.coordinate_offset_y + tree.coordinate_scale
    wx = (x - tree.coordinate_offset_x)/tree.coordinate_scale
    wy = (y - tree.coordinate_offset_y)/tree.coordinate_scale
    vy_left = problem.vy_ref(x_min, y)
    vy_right = problem.vy_ref(x_max, y)
    vy_down = problem.vy_ref(x, y_min)
    vy_up = problem.vy_ref(x, y_max)
    return (1.0-wx)*vy_left + wx*vy_right + (1.0-wy)*vy_down + wy*vy_up
  end
  # Very redundant (just vy<-->vx)
  function vx0(x::Float64, y::Float64)
    x_min = tree.coordinate_offset_x
    y_min = tree.coordinate_offset_y
    x_max = tree.coordinate_offset_x + tree.coordinate_scale
    y_max = tree.coordinate_offset_y + tree.coordinate_scale
    wx = (x - tree.coordinate_offset_x)/tree.coordinate_scale
    wy = (y - tree.coordinate_offset_y)/tree.coordinate_scale
    vx_left = problem.vx_ref(x_min, y)
    vx_right = problem.vx_ref(x_max, y)
    vx_down = problem.vx_ref(x, y_min)
    vx_up = problem.vx_ref(x, y_max)
    return (1.0-wx)*vx_left + wx*vx_right + (1.0-wy)*vx_down + wy*vx_up
  end
  p0 = fill(Float64(0), tree.ne)
  Kcont_dummy = 1.0  # irrelevant since p0 is zero
  v0 = face_field_from_functions(tree, vx0, vy0)
  return vp2sol(tree, v0, p0, Kcont_dummy)
end;

# Overwrite boundary values
function populate_stokes_boundary_velocities!(tree, problem, x)
   for f_id in 1:tree.nf
        if tree.f[f_id].face_type == face_boundary
           if tree.f[f_id].direction == face_vx
             value = problem.vx_ref(get_face_coordinates(tree, f_id)...)
           else
             value = problem.vy_ref(get_face_coordinates(tree, f_id)...)
           end
           x[f_id_to_eqnum(tree, f_id)] = value
        end
   end
   return x
end;
