using Plots
using SimpleTreeMeshes

struct Problem
  name::String
  boundary_conditions::StokesBoundaryConditions
  eta::Function
  eta_min::Float64
  eta_max::Float64
  eta_ref::Float64
  f_x::Function
  f_y::Function
  f_p::Function
  vx_ref::Function
  vy_ref::Function
  p_ref::Function
  coordinate_scale::Float64
  coordinate_offset_x::Float64
  coordinate_offset_y::Float64
end

function example()

    # Define a refinement function
    max_level = 4
    function rf_demo(x, y, d, level)::Bool
        if level >= max_level; return false; end

        # Refine if the 4 corners and center of an element are not all on one side of a circle
        r = 0.25;  cx = 0.7; cy = 0.4
        inside0 = ((x + d/2 - cx)^2 + (y + d/2 - cy)^2 < r^2)
        for (px, py) in ((x, y), (x + d, y), (x, y + d), (x + d, y + d))
            if ((px - cx)^2 + (py - cy)^2 < r^2) != inside0
                return true
            end
        end
        return false
    end;


    # Create the mesh
    mesh = CreateSimpleTreeMesh(rf_demo)



    # Build a Stokes system and solve

    function SolMMSIsoFS()
        # Manufactured Solution with Free-Slip BCs and non-physical / non-zero f_p
        # forcing terms from sympy (see script in this directory)

        eta(x, y)  = 1.0
        eta_min = 1.0
        eta_max = 1.0
        eta_ref = 1.0

        vx_ref(x, y) = sin(pi*x) * cos(pi*y)
        vy_ref(x, y) = sin(2*pi*y) * cos(2*pi*x)
        p_ref(x, y) = 1.0

        f_vx(x, y) = -3.0*pi^2*sin(pi*x)*cos(pi*y) - 4.0*pi^2*sin(2*pi*x)*cos(2*pi*y)
        f_vy(x, y) = -1.0*pi^2*sin(pi*y)*cos(pi*x) - 12.0*pi^2*sin(2*pi*y)*cos(2*pi*x)
        f_p(x, y) =  pi*cos(pi*x)*cos(pi*y) + 2*pi*cos(2*pi*x)*cos(2*pi*y)

        return Problem("MMSIsoFS", stokes_boundary_free_slip, eta, eta_min, eta_max, eta_ref,
                                   f_vx, f_vy, f_p, vx_ref, vy_ref, p_ref, 1.0, 0.0, 0.0)
    end;

    problem = SolMMSIsoFS()

    h_min = 2.0^(-max_level)
    A, b, Kcont = assemble_system(mesh, problem, h_min)
    x = A\b
    v, p = sol2vp(mesh, x, Kcont)

    # Plots
    p_mesh = plot(size=(400,400))
    plot_grid!(p_mesh, mesh)

    p_numbers = plot(size=(400,400))
    plot_grid!(p_numbers, mesh)
    plot_face_numbers!(p_numbers, mesh)
    plot_element_numbers!(p_numbers, mesh, offset=mesh.nf) # offset by the number of faces

    # Plot pressure and velocity
    p_vx = plot(axis=nothing, xaxis=false, yaxis=false, size=(400, 400), aspect_ratio=1)
    plot_face_field!(p_vx, mesh, v, face_vx, minimum(v), maximum(v))
    plot_grid!(p_vx, mesh)
    title!(p_vx, "vx")
    p_vx = add_colorbar(p_vx, minimum(v), maximum(v))

    p_vy = plot(axis=nothing, xaxis=false, yaxis=false, size=(400, 400), aspect_ratio=1)
    plot_face_field!(p_vy, mesh, v, face_vy, minimum(v), maximum(v))
    plot_grid!(p_vy, mesh)
    title!(p_vy, "vy")
    p_vy = add_colorbar(p_vy, minimum(v), maximum(v))

    p_stokes = plot(axis=nothing, xaxis=false, yaxis=false, size=(800, 800), aspect_ratio=1)
    plot_element_field!(p_stokes, mesh, p)
    plot_grid!(p_stokes, mesh, RGB(0.25, 0.25, 0.25))
    plot_averaged_velocity_field!(p_stokes, mesh, v, color=RGB(0.25, 0.25, 0.25))

		l = @layout [a b; _ c; d e]
		p_all = plot(p_mesh, p_numbers, p_vx, p_vy, p_stokes, layout = l, size=(800, 1200))
		display(p_all)

end

example()
