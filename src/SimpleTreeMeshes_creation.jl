### Creation ###################################################################
# p4est is only used in CreateSimpleTreeMesh
using P4est_wrapper

# Some helpers which aren't wrapped yet
const P4EST_MAXLEVEL = Int8(30) # hard-coded default
const P4EST_ROOT_LEN = Int32(1) << P4EST_MAXLEVEL
P4EST_QUADRANT_LEN(level::Int8) = Int32(1) << (P4EST_MAXLEVEL - level) / P4EST_ROOT_LEN;

using MPI
using StaticArrays

function wrap_refinement_function(tree::SimpleTreeMesh, simple_function)
  function wrapped(::Ptr{p4est_t}, which_tree::p4est_topidx_t, quadrant_ptr::Ptr{p4est_quadrant_t})::Cint
        @assert which_tree == 0
        q = unsafe_load(quadrant_ptr)
        x, y = transform_p4est_coordinates(tree, q.x, q.y)
        d = transform_p4est_level_to_element_size(tree, q.level)
        return simple_function(x, y, d, q.level) ? 1 : 0
   end
   return wrapped
end

"""
    Creates a new `SimpleTreMesh` object.

    `refinement_function(x, y, d, level)` should return a boolean value
    indicating if an element with corner at `(x,y)`, size `d`, and
    refinement level `level` should be refined.
"""
function CreateSimpleTreeMesh(refinement_function, extra_uniform_refinement=0;
                              coordinate_scale=1.0, coordinate_offset_x=0.0,
                              coordinate_offset_y=0.0)::SimpleTreeMesh
    # helper functions to decode p4est's compact binary information
    function lnodes_decode(face_code::Int8)::SVector{4, Int8}
        @assert 0 <= face_code <= 0b1111
        if face_code == 0; return SVector{4, Int8}(-1, -1, -1, -1); end
        hang_1 = face_code & 0b0100 != 0
        hang_2 = face_code & 0b1000 != 0
        which_child = face_code & 0b11
        if which_child == 0b00
            return SVector{4, Int8}(hang_1 ? 0 : -1, -1, hang_2 ? 0 : -1, -1)
        elseif which_child == 0b01
            return SVector{4, Int8}(-1, hang_1 ? 0 : -1, hang_2 ? 1 : -1, -1)
        elseif which_child == 0b10
            return SVector{4, Int8}(hang_1 ? 1 : -1, -1, -1, hang_2 ? 0 : -1)
        else
            @assert which_child == 0b11
            return SVector{4, Int8}(-1, hang_1 ? 1 : -1, -1, hang_2 ? 1 : -1)
        end
    end;

    function lnodes_decode_corner(face_code::Int8)::SVector{4, Int8}
        @assert 0 <= face_code <= 0b1111
        if face_code == 0; return SVector{4, Int8}(-1, -1, -1, -1); end
        hang_1 = face_code & 0b0100 != 0
        hang_2 = face_code & 0b1000 != 0
        which_child = face_code & 0b11
        # here we translate hanging faces to hanging corners, using the Z order
        # for the corners
        if which_child == 0b00
            return SVector{4, Int8}(-1, hang_2 ? 1 : -1, hang_1 ? 1 : -1, -1)
        elseif which_child == 0b01
            return SVector{4, Int8}(hang_2 ? 1 : -1, -1, -1, hang_1 ? 1 : -1)
        elseif which_child == 0b10
            return SVector{4, Int8}(hang_1 ? 1 : -1, -1, -1, hang_2 ? 1 : -1)
        else
            @assert which_child == 0b11
            return SVector{4, Int8}(-1, hang_1 ? 1 : -1, hang_2 ? 1 : -1, -1)
        end
    end;


    tree = SimpleTreeMesh()
    tree.coordinate_scale = coordinate_scale
    tree.coordinate_offset_x = coordinate_offset_x
    tree.coordinate_offset_y = coordinate_offset_y

    # wrapped refinement function
    refine_fn_wrapped = wrap_refinement_function(tree, refinement_function)
    refine_fn_c = @cfunction($refine_fn_wrapped, Cint, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant}))

    # Helper for uniform refinement
    function refine_fn_uniform(::Ptr{p4est_t}, ::p4est_topidx_t, ::Ptr{p4est_quadrant_t})::Cint
      return 1
    end
    refine_fn_uniform_c = @cfunction($refine_fn_uniform, Cint, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant}))

    # Initialize p4est, which is only used in this function
    verbosity = P4est_wrapper.SC_LP_ERROR
    #verbosity = P4est_wrapper.SC_LP_DEFAULT # more info, useful for testing
    sc_init(MPI.COMM_WORLD, Cint(true), Cint(true), C_NULL, verbosity)
    p4est_init(C_NULL, verbosity)

    # Create a refined, 2-1 balanced mesh using p4est
    connectivity = p4est_connectivity_new_unitsquare()
    forest = p4est_new(MPI.COMM_WORLD, connectivity, 0, C_NULL, C_NULL);
    p4est_refine(forest, 1, refine_fn_c, C_NULL)  # second argument is recursive refinement
    p4est_balance(forest, P4est_wrapper.P4EST_CONNECT_FULL, C_NULL)
    for  _ in 1:extra_uniform_refinement
      p4est_refine(forest, 0, refine_fn_uniform_c, C_NULL)  # second argument is recursive refinement
    end

    ghost = p4est_ghost_new(forest, P4est_wrapper.P4EST_CONNECT_FULL)

    # For convenience, use p4est's numbering routines, even though they
    # are specialized for finite element methods.
    lnodes_face_ptr = p4est_lnodes_new(forest, ghost, -1); # one node per non-hanging face
    lnodes_face = unsafe_load(lnodes_face_ptr)

    # Assuming a single tree on one rank, extract the number of local elements
    # Note that we do not do this for face and corner counts, as the
    # "lnodes" functions we're using as helpers do not include hanging faces and corners
    # (For maximum performance, one might want to write custom versions of these
    # functions which exactly match our use case, where these hanging entities are not
    # always eliminated, as they are in the finite element setting which motivates p4est)
    tree.ne = lnodes_face.num_local_elements

    # Iterate over all elements with p4est_iterate to populate e
    tree.e = Vector{ElementInfo}(undef, tree.ne)
    function volume_callback_e(info_ptr::Ptr{p4est_iter_volume_info_t}, ::Ptr{Cvoid})::Cvoid
        info = unsafe_load(info_ptr);
        q = unsafe_load(info.quad)
        i = info.quadid + 1  # p4est is 0-based, we are 1-based
        tree.e[i] = ElementInfo(q.x, q.y, q.level)
        return nothing
    end
    volume_callback_c = @cfunction($volume_callback_e, Cvoid, (Ptr{p4est_iter_volume_info_t}, Ptr{Cvoid}));
    p4est_iterate(forest, C_NULL, C_NULL, volume_callback_c, C_NULL, C_NULL)

    # Iterate again (no need for p4est_iterate) and collect how many face dofs (1 or 3)
    # correspond to each non-hanging face.
    nf_nohang = lnodes_face.num_local_nodes
    f_per_f_nohang = fill(Int8(1), nf_nohang)
    nf_count = nf_nohang
    for i = 1:tree.ne
        code = unsafe_load(lnodes_face.face_code, i)
        decoded = lnodes_decode(code)
        for v = 1:lnodes_face.vnodes
            f_nohang = unsafe_load(lnodes_face.element_nodes, (lnodes_face.vnodes * (i-1)) + v) + 1
            if decoded[v] >= 0
                if f_per_f_nohang[f_nohang] == 1
                  f_per_f_nohang[f_nohang] = 3
                  nf_count += 2
                end
            end
        end
    end
    tree.nf = nf_count

    # condense f_per_nohang to get f_nohang_to_f map
    f_nohang_to_f = Vector{Int32}(undef, nf_nohang)
    j = 1
    for i = 1:nf_nohang
        f_nohang_to_f[i] = j
        j += f_per_f_nohang[i]
    end

    # Iterate over all elements again
    tree.f = Vector{FaceInfo}(undef, tree.nf)
    tree.e2f = fill(-1, (FACES_PER_ELEMENT, tree.ne))
    tree.f2e = fill(-1, (ELEMENTS_PER_FACE, tree.nf))
    for i = 1:tree.ne
        code = unsafe_load(lnodes_face.face_code, i)
        decoded = lnodes_decode(code)
        for v = 1:lnodes_face.vnodes
            f_nohang = unsafe_load(lnodes_face.element_nodes, (lnodes_face.vnodes * (i-1)) + v) + 1
            off = 1 + decoded[v]
            if decoded[v] >= 0
                if off == 1; face_type = face_first; end
                if off == 2; face_type = face_second; end
            else
              face_type = f_per_f_nohang[f_nohang] == 1 ? face_regular : face_big
            end
            f_id = f_nohang_to_f[f_nohang] + off
            tree.e2f[v, i] = f_id
            direction = v in [LEFT, RIGHT] ? face_vx : face_vy
            if v in [RIGHT, UP]
              tree.f2e[1, f_id] = i
            else
              @assert v in [LEFT, DOWN]
              tree.f2e[2, f_id] = i
            end
            tree.f[f_id] = FaceInfo(face_type, direction)
        end
    end

    # Iterate again to determine boundary faces with populated f2e.
    # This is ugly as it involves overwriting everything
    for f_id = 1:tree.nf
      if tree.f[f_id].face_type == face_regular && -1 in tree.f2e[:, f_id]
        tree.f[f_id] = FaceInfo(face_boundary, tree.f[f_id].direction)
      end
    end

    # Repeat a similar procedure for the corners
    lnodes_corner_ptr = p4est_lnodes_new(forest, ghost, 1); # one node per non-hanging corner
    lnodes_corner = unsafe_load(lnodes_corner_ptr)

    # Iterate once again to populate a map from non-hanging corners to elements
    nc_nohang = lnodes_corner.num_local_nodes
    c_nohang2e = fill(-1, (ELEMENTS_PER_CORNER, nc_nohang))
    for i = 1:tree.ne
        code = unsafe_load(lnodes_corner.face_code, i)
        decoded_corner = lnodes_decode_corner(code)

        # Note that these will be from the *parent* for hanging nodes.
        for v in [DOWN_LEFT, DOWN_RIGHT, UP_LEFT, UP_RIGHT]
            if decoded_corner[v] == -1
                c_nohang_id = unsafe_load(lnodes_corner.element_nodes, (4 * (i-1)) + v) + 1
                c_nohang2e[opp_diag(v), c_nohang_id] = i
            end
        end
    end

    # Now, we can iterate through the nonhanging corners and determine
    # which have neighboring hanging corners
    # Populate c_nohang_to_c
    # Partially populate c2h FIXME: not sure we actually need all this
    c_nohang2h = fill(-1, (4, nc_nohang)) # likely inefficient. left-right-down-up convention
    c_nohang_to_c = fill(-1, nc_nohang)
    idx = 0
    for i = 1:nc_nohang
        idx += 1
        c_nohang_to_c[i] = idx

        # We are interested in checking (only) to see if this non-hanging node
        # has a hanging node above or to the right. So we check for level differences
        # in the neighboring elements in those directions.
        # Note: it might actually be better to do down and left, so that hanging nodes
        #       would be numbered after both neighboring non-hanging nodes
        e_id_up_right = c_nohang2e[UP_RIGHT, i]
        e_id_down_right = c_nohang2e[DOWN_RIGHT, i]
        e_id_up_left = c_nohang2e[UP_LEFT, i]
        if e_id_up_right != -1
            if e_id_up_left != -1
                if tree.e[e_id_up_right].level != tree.e[e_id_up_left].level;
                    idx += 1;
                    c_nohang2h[UP, i] = idx;
                end
            end
            if e_id_down_right != -1
                if tree.e[e_id_up_right].level != tree.e[e_id_down_right].level
                    idx += 1
                    c_nohang2h[RIGHT, i] = idx;
                end
            end
        end
    end
    tree.nc = idx

    # Prepare to populate corner types by flagging hanging corners
    # FIXME: this seems a bit hacky, still
    tree.c = fill(CornerInfo(corner_undefined, 0), tree.nc)
    for i=1:nc_nohang
      for v in [UP, RIGHT]
        idx = c_nohang2h[v ,i]
        if idx != -1
          tree.c[idx] = CornerInfo(corner_hanging, 0)
        end
      end
    end

    # Iterate over all elements to populate e2c and c2e map
    tree.e2c = fill(-1, (4, tree.ne))
    tree.c2e = fill(-1, (ELEMENTS_PER_CORNER, tree.nc)) # For hanging corners, two on the "big" side are -1
    for i = 1:tree.ne
        code = unsafe_load(lnodes_corner.face_code, i)
        decoded_corner = lnodes_decode_corner(code)
        which_child = code & 0b11

        # Note that these will be from the *parent* for hanging nodes.
        c_nohang_id = [-1, -1, -1, -1]
        for v = 1:4
            c_nohang_id[v] =  unsafe_load(lnodes_corner.element_nodes, (4 * (i-1)) + v) + 1
        end

        for v in [DOWN_LEFT, DOWN_RIGHT, UP_LEFT, UP_RIGHT]

            # If hanging, check each of the pair of c_nohang_id's in c2h.
            # if not found, add to the second of the pair.
            # Note that we assume we have already populated all of the "right" and "up" entries
            if decoded_corner[v] == -1
                c_id_full = c_nohang_to_c[c_nohang_id[v]]
                tree.e2c[v, i] = c_id_full
                tree.c2e[opp_diag(v), c_id_full] = i
            else
                if which_child == 0b00 # downleft child
                    if v == UP_LEFT
                        v1 = DOWN_LEFT; dir1 = UP;
                        v2 = UP_LEFT; dir2 = DOWN;
                    else
                        @assert v == DOWN_RIGHT # only two possible hangers
                        v1 = DOWN_LEFT; dir1 = RIGHT;
                        v2 = DOWN_RIGHT; dir2 = LEFT;
                    end
                elseif which_child == 0b01 # DOWNRIGHT child
                    if v == DOWN_LEFT
                        v1 = DOWN_LEFT; dir1 = RIGHT;
                        v2 = DOWN_RIGHT; dir2 = LEFT;
                    else
                        @assert v == UP_RIGHT
                        v1 = DOWN_RIGHT; dir1 = UP;
                        v2 = UP_RIGHT; dir2 = DOWN;
                    end
                elseif which_child == 0b10 # UPLEFT child
                    if v == DOWN_LEFT
                        v1 = DOWN_LEFT; dir1 = UP;
                        v2 = UP_LEFT; dir2 = DOWN;
                    else
                        @assert v == UP_RIGHT
                        v1 = UP_LEFT; dir1 = RIGHT;
                        v2 = UP_RIGHT; dir2 = LEFT;
                    end
                else # UPRIGHT child
                    @assert which_child == 0b11
                    if v == UP_LEFT
                        v1 = UP_LEFT; dir1 = RIGHT;
                        v2 = UP_RIGHT; dir2 = LEFT;
                    else
                        @assert v == DOWN_RIGHT
                        v1 = DOWN_RIGHT; dir1 = UP;
                        v2 = UP_RIGHT; dir2 = DOWN;
                    end
                end
                check1 = c_nohang2h[dir1, c_nohang_id[v1]]
                if check1 != -1
                    c_id_hanging = check1

                    # Fill in a missing entry of c_nohang2h
                    c_nohang2h[dir2, c_nohang_id[v2]] = c_id_hanging
                else
                    # All hanging nodes should have already been accounted for
                    check2 = c_nohang2h[dir2, c_nohang_id[v2]]
                    @assert check2 != -1
                    c_id_hanging = check2
                end
                tree.e2c[v, i] = c_id_hanging
                tree.c2e[opp_diag(v), c_id_hanging] = i
            end
        end
    end

    p4est_lnodes_destroy(lnodes_corner_ptr)
    p4est_lnodes_destroy(lnodes_face_ptr)
    p4est_ghost_destroy(ghost)
    p4est_destroy(forest)
    p4est_connectivity_destroy(connectivity)
    sc_finalize()
    if (MPI.Initialized() && !isinteractive())
        MPI.Finalize()
    end

    # Complete corner types
    for i = 1:tree.nc
      is_hanging = tree.c[i].corner_type == corner_hanging
      is_boundary = !is_hanging && -1 in tree.c2e[:, i]

      if is_boundary
        c2e_local = tree.c2e[:, i]
        present = findall(!=(-1), c2e_local)
        if present == [1]
            info = CORNER_BOUNDARY_UP_RIGHT
        elseif present == [2]
            info = CORNER_BOUNDARY_UP_LEFT
        elseif present == [3]
            info = CORNER_BOUNDARY_DOWN_RIGHT
        elseif present == [4]
            info = CORNER_BOUNDARY_DOWN_LEFT
        else
            @assert length(present) == 2
            level1 = tree.e[c2e_local[present[1]]].level
            level2 = tree.e[c2e_local[present[2]]].level
            if level1 < level2
                level_bits = 0b0100
            elseif level1 > level2
                level_bits = 0b1000
            else
                level_bits = 0b0000
            end
            if present == [2, 4] # left boundary
                info = level_bits + 0b00
            elseif present == [1, 3] # right boundary
                info = level_bits + 0b01
            elseif present == [3, 4] # bottom boundary
                info = level_bits + 0b10
            elseif present == [1, 2] # top boundary
                info = level_bits + 0b11
            else
                @assert false
            end
        end
      else
        levels = [-1, -1, -1, -1]
        min_level = P4EST_MAXLEVEL
        for v = 1:4
          e_id = tree.c2e[v, i]
          if e_id != -1
            level = tree.e[e_id].level
            min_level = min(min_level, level)
            levels[v] = level
          end
        end

        if is_hanging
          min_level -= 1 # count the level of the unreferenced larger neighbor
          @assert min_level >= 0
        end

        for v = 1:4
          if levels[v] == -1
            levels[v]= 0
          else
            levels[v] -= min_level
          end
        end
        info = 0b1000 * levels[1] + 0b0100 * levels[2] + 0b0010 * levels[3] + 0b0001 * levels[4]
      end

      if is_hanging
        tree.c[i] = CornerInfo(corner_hanging, info)
      elseif is_boundary
        tree.c[i] = CornerInfo(corner_boundary, info)
      else
        tree.c[i] = CornerInfo(corner_regular, info)
      end
    end

    # Populate corner-face maps.
    # This is done after all p4est/MPI operations, so could be done
    # with a separate function, lazily.
    tree.f2c = fill(-1, (CORNERS_PER_FACE, tree.nf))
    tree.c2f = fill(-1, (FACES_PER_CORNER, tree.nc))
    for f_id in 1:tree.nf
      if tree.f[f_id].direction == face_vx
        e_id_left = tree.f2e[1, f_id]
        if  e_id_left != -1 # left element exists
          c_id_down = tree.e2c[DOWN_RIGHT, e_id_left]
          c_id_up = tree.e2c[UP_RIGHT, e_id_left]
        else
          @assert tree.f2e[2, f_id] != -1 # right element exists
          e_id_right = tree.f2e[2, f_id]
          c_id_down = tree.e2c[DOWN_LEFT, e_id_right]
          c_id_up = tree.e2c[UP_LEFT, e_id_right]
        end
        tree.f2c[1, f_id] = c_id_down
        tree.f2c[2, f_id] = c_id_up
        # c2f prefers "face_big" over "face_first" and "face_second"
        if tree.f[f_id].face_type ∉ (face_first, face_second) || tree.c2f[DOWN, c_id_up] == -1
          tree.c2f[DOWN, c_id_up] = f_id
        end
        if tree.f[f_id].face_type ∉ (face_first, face_second) || tree.c2f[UP, c_id_down] == -1
          tree.c2f[UP, c_id_down] = f_id
        end
      else
        @assert tree.f[f_id].direction == face_vy
        e_id_down = tree.f2e[1, f_id]
        if e_id_down != -1 # down element exists
          c_id_left = tree.e2c[UP_LEFT, e_id_down]
          c_id_right = tree.e2c[UP_RIGHT, e_id_down]
        else
          @assert tree.f2e[2, f_id] != -1 # up element exists
          e_id_up = tree.f2e[2, f_id]
          c_id_left = tree.e2c[DOWN_LEFT, e_id_up]
          c_id_right = tree.e2c[DOWN_RIGHT, e_id_up]
        end
        tree.f2c[1, f_id] = c_id_left
        tree.f2c[2, f_id] = c_id_right
        # c2f prefers "face_big" over "face_first" and "face_second"
        if tree.f[f_id].face_type ∉ (face_first, face_second) || tree.c2f[LEFT, c_id_right] == -1
          tree.c2f[LEFT, c_id_right] = f_id
        end
        if tree.f[f_id].face_type ∉ (face_first, face_second) || tree.c2f[RIGHT, c_id_left] == -1
          tree.c2f[RIGHT, c_id_left] = f_id
        end
      end
    end

    return tree
end;
