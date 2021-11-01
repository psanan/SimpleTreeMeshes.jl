### Embedding/Coordinates ######################################################

function transform_p4est_level_to_element_size(tree::SimpleTreeMesh, level)::Float64
    return tree.coordinate_scale * Float64(P4EST_QUADRANT_LEN(level))
end

function get_element_size(tree::SimpleTreeMesh, e_id::Int)
    return transform_p4est_level_to_element_size(tree, tree.e[e_id].level)
end

function get_face_size(tree::SimpleTreeMesh, f_id::Int)
  e_id_1 = tree.f2e[1, f_id]
  if e_id_1 != -1
    return get_element_size(tree, e_id_1)
  else
    return get_element_size(tree, tree.f2e[2, f_id])
  end
end

function get_corner_coordinates(tree::SimpleTreeMesh, c_id)
    xc = yc = 0
    for v = 1:CORNERS_PER_ELEMENT
        e_id = tree.c2e[v, c_id]
        if e_id != -1
            x, y = get_element_corner_coordinates(tree, e_id)
            d = get_element_size(tree, e_id)
            if v == DOWN_LEFT
                xc = x+d; yc = y+d;
            elseif v == DOWN_RIGHT
                xc = x; yc = y+d;
            elseif v == UP_LEFT
                xc = x+d; yc = y;
            else
                @assert v == UP_RIGHT
                xc = x; yc = y;
            end
            break
        end
    end
    return xc,yc
end

function get_face_coordinates(tree::SimpleTreeMesh, f_id)
    direction = tree.f[f_id].direction
    if direction == face_vx
        e_id_left = tree.f2e[1, f_id]
        if e_id_left != -1
            e_id = e_id_left
            which_face = RIGHT
        else
            @assert tree.f2e[2, f_id] != -1
            e_id = tree.f2e[2, f_id]
            which_face = LEFT
        end
    else
        @assert direction == face_vy
        e_id_down = tree.f2e[1, f_id]
        if e_id_down != -1
            e_id = e_id_down
            which_face = UP
        else
            @assert tree.f2e[2, f_id] != -1
            e_id = tree.f2e[2, f_id]
            which_face = DOWN
        end
    end
    x_e, y_e = get_element_corner_coordinates(tree, e_id)
    d = get_element_size(tree, e_id)
    if which_face == LEFT
        x = x_e
        y = y_e + d/2
    elseif which_face == RIGHT
        x = x_e + d
        y = y_e + d/2
    elseif which_face == DOWN
        x = x_e + d/2
        y = y_e
    else
        @assert which_face == UP
        x = x_e + d/2
        y = y_e + d
    end
    return x, y
end

function transform_p4est_coordinates(tree::SimpleTreeMesh, x_p4est, y_p4est)::NTuple{2,Float64}
    x = tree.coordinate_scale * Float64(x_p4est) / P4EST_ROOT_LEN + tree.coordinate_offset_x
    y = tree.coordinate_scale * Float64(y_p4est) / P4EST_ROOT_LEN + tree.coordinate_offset_y
    return x, y
end

function get_element_corner_coordinates(tree::SimpleTreeMesh, e_id)
    return transform_p4est_coordinates(tree, tree.e[e_id].x, tree.e[e_id].y)
end

function get_element_coordinates(tree::SimpleTreeMesh, e_id)
    x, y = get_element_corner_coordinates(tree, e_id)
    d = get_element_size(tree, e_id)
    return x+d/2, y+d/2
end

### Topology Helpers ###########################################################

function hanging_corner_get_big_element(tree::SimpleTreeMesh, c_id)
    @assert tree.c[c_id].corner_type == corner_hanging
    info = tree.c[c_id].info
    if info == CORNER_HANG_LEFT
        # Note: this case and all the others could be much more concise without the asserts
        e_id_up_right = tree.c2e[UP_RIGHT, c_id]
        @assert e_id_up_right != -1
        c_id_up = tree.e2c[UP_LEFT, e_id_up_right]
        e_id_left = tree.c2e[DOWN_LEFT, c_id_up]
        @assert e_id_left != -1
        return e_id_left
    elseif info == CORNER_HANG_RIGHT
        e_id_up_left = tree.c2e[UP_LEFT, c_id]
        @assert e_id_up_left != -1
        c_id_up = tree.e2c[UP_RIGHT, e_id_up_left]
        e_id_right = tree.c2e[DOWN_RIGHT, c_id_up]
        @assert e_id_right != -1
        return e_id_right
    elseif info == CORNER_HANG_DOWN
        e_id_up_right = tree.c2e[UP_RIGHT, c_id]
        @assert e_id_up_right != -1
        c_id_right = tree.e2c[DOWN_RIGHT, e_id_up_right]
        e_down = tree.c2e[DOWN_LEFT, c_id_right]
        @assert e_down != -1
        return e_down
    else
        @assert info == CORNER_HANG_UP
        e_id_down_right = tree.c2e[DOWN_RIGHT, c_id]
        @assert e_id_down_right != 0
        c_id_right = tree.e2c[UP_RIGHT, e_id_down_right]
        e_id_up = tree.c2e[UP_LEFT, c_id_right]
        @assert e_id_up != -1
        return e_id_up
    end
end

function c2f_small(tree::SimpleTreeMesh, loc, c_id::Int)
    f_id = tree.c2f[loc, c_id]
    if tree.f[f_id].face_type == face_big
        if loc in (RIGHT, UP)
            f_id += 1
            @assert tree.f[f_id].face_type == face_first
        else
            @assert loc in (LEFT, DOWN)
            f_id += 2
            @assert tree.f[f_id].face_type == face_second
        end
    end
    return f_id
end
