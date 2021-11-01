# Note: we are probably reinventing wheels here. Consider using something like Makie,
#       e.g. as is done at https://github.com/tduretz/FCFV_Geodynamics.jl/blob/main/VisuFCFV.jl

### Plotting ###################################################################
# We define some helpers to plot our mesh, numberings on it, and fields living on elements, faces, and corners.
#For fields, we aim not to plot nice-looking interpolations of the fields, but rather collections of cells corresponding
#to the control volumes associated with degrees of freedom in the discretization of the Stokes equations themselves.
#This should be suitable to most directly show the performance of the method. Note that here we have chosen volumes which were easy to compute, not necessarily those which are used in particular methods.iWe define some helpers to plot our mesh, numberings on it, and fields living on elements, faces, and corners.
#For fields, we aim not to plot nice-looking interpolations of the fields, but rather collections of cells corresponding
#to the control volumes associated with degrees of freedom in the discretization of the Stokes equations themselves.
#This should be suitable to most directly show the performance of the method. Note that here we have chosen volumes which were easy to compute, not necessarily those which are used in particular methods.

using Plots

struct plot_grid_data; rect_x; rect_y; end

function plot_grid!(p, tree::SimpleTreeMesh, color=:black)
    data = plot_grid_data(Vector{Float16}(), Vector{Float16}())
    for e_id = 1:tree.ne
        x, y = get_element_corner_coordinates(tree, e_id)
        d = get_element_size(tree, e_id)
        append!(data.rect_x, [x, x+d, x+d, x, x, NaN])
        append!(data.rect_y, [y, y, y+d, y+d, y, NaN])
    end
    plot!(p, data.rect_x, data.rect_y, color=color, label=nothing)
    return nothing
end;

function plot_text!(p, x, y, value, level, color=:blue)
    text_size = max(4,Int(floor(16 - 2.7*level)))
    annotate!(p, x, y, text(value, color, text_size))
end;

const PLOT_TEXT_MAX_LEVEL = 4;

function plot_element_numbers!(p, tree::SimpleTreeMesh; offset = 0, max_level = PLOT_TEXT_MAX_LEVEL)
   for e_id = 1:tree.ne
        level = tree.e[e_id].level
        if level <= max_level
            x, y = get_element_corner_coordinates(tree, e_id)
            d = get_element_size(tree, e_id)
            plot_text!(p, x + d/2, y + d/2, e_id + offset, level, :red)
        end
    end
end;

# FIXME: now that we have a function to get coords for faces, and more face info, this should iterate over faces, plotting at the center, except for split faces where we do the offsets as below:
function plot_face_numbers!(p, tree::SimpleTreeMesh; offset = 0, max_level = PLOT_TEXT_MAX_LEVEL)
    for e_id = 1:tree.ne
        for v = 1:FACES_PER_ELEMENT
            level = tree.e[e_id].level
            if level <= max_level
                x, y = get_element_corner_coordinates(tree, e_id)
                d = get_element_size(tree, e_id)
                if v == LEFT
                    x_plot, y_plot = x + d/8, y + d/2
                elseif v == RIGHT
                    x_plot, y_plot = x + d*7/8, y + d/2
                elseif v == DOWN
                    x_plot, y_plot = x + d/2, y + d/8
                else
                    @assert v == UP
                    x_plot, y_plot = x + d/2, y + d*7/8
                end
                value = tree.e2f[v, e_id] + offset
                plot_text!(p, x_plot, y_plot, value, level, :green)
            end
        end
    end
end;

function plot_corner_numbers!(p, tree::SimpleTreeMesh; offset = 0, max_level = PLOT_TEXT_MAX_LEVEL)
    for i = 1:tree.nc
        min_level = P4EST_MAXLEVEL
        for v=1:CORNERS_PER_ELEMENT
            e_id = tree.c2e[v, i]
            if e_id != -1
                level = tree.e[e_id].level
                min_level = minimum([level, min_level])
            end
        end
        if min_level <= max_level
            x_plt,y_plt = get_corner_coordinates(tree, i)
            type = tree.c[i].corner_type
            if type == corner_regular
               color = :blue
            elseif type == corner_boundary
               color = :magenta
            elseif type == corner_undefined
               color = :red
            else
              @assert type == corner_hanging
              color = :blueviolet
            end
          plot_text!(p, x_plt, y_plt, i + offset, min_level, color)
        end
    end
end;

function colormap(value)
  if !(0 <= value <= 1)
    println("ERROR value ", value, " is not in [0,1]")
  end
  @assert 0 <= value <= 1
  return RGB(value, 100/255, 100/255)
end;

# FIXME probably doesn't look right when not using default coords
function colorbar!(p, min_value=0.0, max_value=1.0; n_segments=30)
  if min_value == max_value; max_value = min_value + 1; end
  @assert min_value < max_value
  h = (max_value - min_value) / n_segments
  for i = 1:n_segments
    c = colormap((i-1)/(n_segments -1))
    plot!(
          p,
          Shape([0, 1,  1, 0], min_value .+ (i-1)*h .+ [0, 0, h, h]),
          fill = c,
          label = nothing,
          linewidth = 1,
          linecolor = c,
         )
  end
  if min_value < 0 < max_value
     ticks = [min_value, 0, max_value]
  else
     ticks = [min_value, max_value]
  end
  plot!(p, yticks=round.(ticks, sigdigits=3))
end;

function add_colorbar(p, min_value, max_value)
    plot_size_x, plot_size_y = p[:size]
    colorbar_width = 20
    p_colorbar = plot(axis=nothing, xaxis=false, yaxis=false, size=(colorbar_width, plot_size_y))
    colorbar!(p_colorbar, min_value, max_value)
    layout = @layout [a{0.05w} b]
    combined_size = (colorbar_width + plot_size_x, plot_size_y)
    p_combined = plot(p_colorbar, p, layout=layout, size=combined_size)
    return p_combined
end;

function plot_rect!(p, x, y, w, h, value, min_value, max_value)
    value_scaled = (value - min_value)/(max_value - min_value)
    c = colormap(value_scaled)
    plot!(
        p,
        Shape(x .+ [0,w,w,0], y .+ [0,0,h,h]),
        fill=c,
        linecolor=c,
        linewidth=1, # FIXME: this can cause problems!
        label=nothing
    )
end;

function plot_element_field!(p, tree::SimpleTreeMesh, field, min_value=nothing, max_value=nothing)
  @assert length(field) == tree.ne
  if min_value == nothing; min_value = minimum(field); end
  if max_value == nothing; max_value = maximum(field); end
  if min_value == max_value; max_value = min_value + 1; end
  for e_id = 1:tree.ne
        x, y = get_element_corner_coordinates(tree, e_id)
        d = get_element_size(tree, e_id)
        plot_rect!(p, x, y, d, d, field[e_id], min_value, max_value)
    end
end;

# Note: this currently makes the (easiest-to-implement) choice
# that "big" edges have their own data, and we plot that data as with all
# other values, as half an element. However, one might also want the option
# to plot nothing for these edges and instead divide the corresponding
# region between the two "small" edges, which corresponds to a finite volume
# method.
function plot_face_field!(p, tree::SimpleTreeMesh, field, direction::FaceDirection=face_vx, min_value=nothing, max_value=nothing)
    @assert direction == face_vx || direction == face_vy
    @assert length(field) == tree.nf
    if min_value == nothing; min_value = minimum(field); end
    if max_value == nothing; max_value = maximum(field); end
    if min_value == max_value; max_value = min_value + 1; end
    for e_id = 1:tree.ne
        x, y = get_element_corner_coordinates(tree, e_id)
        d = get_element_size(tree, e_id)

        faces = direction == face_vx ? [LEFT, RIGHT] : [DOWN, UP]
        for v in faces
            f_id = tree.e2f[v, e_id]
            value = field[f_id]
            if v == LEFT
              plot_rect!(p, x, y, d/2, d, value, min_value, max_value)
            elseif v == RIGHT
              plot_rect!(p, x+d/2, y, d/2, d, value, min_value, max_value)
            elseif v == DOWN
              plot_rect!(p, x, y, d, d/2, value, min_value, max_value)
            else
                @assert v == UP
                plot_rect!(p, x, y+d/2, d, d/2, value, min_value, max_value)
            end
        end
    end
end;

function plot_corner_field!(p, tree::SimpleTreeMesh, field, min_value=nothing, max_value=nothing)
    @assert length(field) == tree.nc
    if min_value == nothing; min_value = minimum(field); end
    if max_value == nothing; max_value = maximum(field); end
    if min_value == max_value; max_value = min_value + 1; end
    # Perform two passes, saving big element next to hanging corners for last
    for c_id = 1:tree.nc
        for v in [DOWN_LEFT, DOWN_RIGHT, UP_LEFT, UP_RIGHT]
            e_id = tree.c2e[v, c_id]
            if e_id != -1
                x, y = get_element_corner_coordinates(tree, e_id)
                d = get_element_size(tree, e_id)
                if v == DOWN_LEFT
                    x_plt, y_plt = x + d/2, y + d/2
                elseif v == DOWN_RIGHT
                    x_plt, y_plt = x, y + d/2
                elseif v == UP_LEFT
                    x_plt, y_plt = x + d/2, y
                else
                    @assert v == UP_RIGHT
                    x_plt, y_plt = x, y
                end
                plot_rect!(p, x_plt, y_plt, d/2, d/2, field[c_id], min_value, max_value)
            end
        end
    end

    for c_id = 1:tree.nc
        if tree.c[c_id].corner_type == corner_hanging
            e_id = hanging_corner_get_big_element(tree, c_id)
            x, y = get_element_corner_coordinates(tree, e_id)
            d = get_element_size(tree, e_id)
            info = tree.c[c_id].info
            if info == CORNER_HANG_LEFT
                x_plt, y_plt = x + d/2, y + d/4
            elseif info == CORNER_HANG_RIGHT
                x_plt, y_plt = x, y + d/4
            elseif info == CORNER_HANG_DOWN
                x_plt, y_plt = x + d/4, y + d/2
            else
                @assert info == CORNER_HANG_UP
                x_plt, y_plt = x + d/4, y
            end
            plot_rect!(p, x_plt, y_plt, d/2, d/2, field[c_id], min_value, max_value)
        end
    end

end

function arrow0!(p, x, y, u, v; as=0.2, lc=:black)
    nuv = sqrt(u^2 + v^2)
    v1, v2 = [u;v] / nuv,  [-v;u] / nuv
    v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
    v5 = v4 - 2*(v4'*v2)*v2
    v4, v5 = as*nuv*v4, as*nuv*v5
    plot!(p, [x,x+u], [y,y+v], lc=lc, label=nothing)
    plot!(p, [x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lc=lc, label=nothing)
    plot!(p, [x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lc=lc, label=nothing)
end

# Plot directions of averaged velocity field (length based on cell size, not
# velocity magnitude)
function plot_averaged_velocity_field!(p, tree::SimpleTreeMesh, field; color=:black)
    @assert length(field) == tree.nf
    for e_id = 1:tree.ne
        xc, yc = get_element_coordinates(tree, e_id)
        d = get_element_size(tree, e_id)

        vx = (field[tree.e2f[LEFT, e_id]] + field[tree.e2f[RIGHT, e_id]])/2
        vy = (field[tree.e2f[DOWN, e_id]] + field[tree.e2f[UP, e_id]])/2

        scaling = 0.7 * (d/2) / sqrt(vx^2 + vy^2)
        arrow_x = [xc, xc + scaling * vx]
        arrow_y = [yc, yc + scaling * vy]
        arrowsize = 0.00001
        arrow0!(p, xc, yc, scaling*vx, scaling*vy, lc=color)
    end
end;

