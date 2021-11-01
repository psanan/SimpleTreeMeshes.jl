function element_field_from_function(tree::SimpleTreeMesh, point_function)
    field = Vector{Float64}(undef, tree.ne)
    for i = 1:tree.ne
        x = tree.coordinate_scale * Float64(tree.e[i].x) / P4EST_ROOT_LEN + tree.coordinate_offset_x
        y = tree.coordinate_scale * Float64(tree.e[i].y) / P4EST_ROOT_LEN + tree.coordinate_offset_y
        d = tree.coordinate_scale * Float64(P4EST_QUADRANT_LEN(tree.e[i].level))
        field[i] = point_function(x + d/2, y + d/2)
    end
    return field
end;

function face_field_from_functions(tree::SimpleTreeMesh, point_function_vx, point_function_vy)
    field = Vector{Float64}(undef, tree.nf)
    # Note that this redundantly computes on shared edges.
    for i = 1:tree.ne
        x = tree.coordinate_scale * Float64(tree.e[i].x) / P4EST_ROOT_LEN + tree.coordinate_offset_x
        y = tree.coordinate_scale * Float64(tree.e[i].y) / P4EST_ROOT_LEN + tree.coordinate_offset_y
        d = tree.coordinate_scale * Float64(P4EST_QUADRANT_LEN(tree.e[i].level))

        for v in [LEFT, RIGHT]
            f_id = tree.e2f[v, i]
            if v == LEFT
                x_pt, y_pt = x, y+d/2
            else
              @assert v == RIGHT
              x_pt, y_pt = x+d, y+d/2
            end
            field[f_id] = point_function_vx(x_pt, y_pt)
        end
        for v in [DOWN, UP]
            f_id = tree.e2f[v, i]
            if v == DOWN
                x_pt, y_pt = x+d/2, y
            else
                @assert v == UP
                x_pt, y_pt, = x+d/2, y+d
            end
            field[f_id] = point_function_vy(x_pt, y_pt)
        end
    end
    return field
end;

function corner_field_from_function(tree::SimpleTreeMesh, point_function)
    field = Vector{Float64}(undef, tree.nc)
    for c_id = 1:tree.nc
        x, y = get_corner_coordinates(tree, c_id)
        field[c_id] = point_function(x, y)
    end
    return field
end;

function L1_norm_face(tree::SimpleTreeMesh, v)
    @assert length(v) == tree.nf

    # FIXME - naive for now
    norm = 0.0
    for f_id = 1:tree.nf
        da = 1.0 # FIXME placeholder
        norm += da * abs(v[f_id])
    end

    return norm
end;

function L1_norm_element(tree::SimpleTreeMesh, x)
    @assert length(x) == tree.ne

    norm = 0.0
    a_total = 0.0 # just for assertion

    # 1-point quadrature
    for i = 1:tree.ne
        d = get_element_size(tree, i)
        da = d*d
        norm += da * abs(x[i])
        a_total += da # just for assertion
    end
    @assert abs(a_total - (tree.coordinate_scale^2)) < 1e-8

    return norm
end;

function L2_norm_corner(tree::SimpleTreeMesh, x)
    @assert length(x) == tree.nc

    # FIXME - naive for now
    norm_squared = 0.0
    for c_id = 1:tree.nc
        da = 1.0/tree.nc # FIXME placeholder
        norm_squared += da * x[c_id]^2
    end

    return sqrt(norm_squared)
end;

function L2_norm_face(tree::SimpleTreeMesh, v)
    @assert length(v) == tree.nf

    # FIXME - naive for now
    norm_squared = 0.0
    for f_id = 1:tree.nf
        da = 1.0 # FIXME placeholder
        norm_squared += da * v[f_id]^2
    end

    return sqrt(norm_squared)
end;

function L2_norm_element(tree::SimpleTreeMesh, x)
    @assert length(x) == tree.ne

    norm_squared = 0.0
    a_total = 0.0 # just for assertion

    # 1-point quadrature
    for i = 1:tree.ne
        d = get_element_size(tree, i)
        da = d*d
        norm_squared += da * x[i]^2
        a_total += da # just for assertion
    end
    @assert abs(a_total - (tree.coordinate_scale^2)) < 1e-8

    return sqrt(norm_squared)
end;

function face_field_average_to_element_field(tree::SimpleTreeMesh, face_field, direction::FaceDirection)
    @assert length(face_field) == tree.nf
    element_field = Vector{Float64}(undef, tree.ne)
    dir1 = direction==face_vx ? LEFT : DOWN
    dir2 = direction==face_vx ? RIGHT : UP
    for e_id = 1:tree.ne
        f_id1 = tree.e2f[dir1, e_id]
        f_id2 = tree.e2f[dir2, e_id]
        element_field[e_id] = 0.5 * (face_field[f_id1] + face_field[f_id2])
   end
   return element_field
end;
