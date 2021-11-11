module SimpleTreeMeshes

export SimpleTreeMesh, CreateSimpleTreeMesh
export plot_grid_data, plot_grid!, plot_element_numbers!, plot_face_numbers!, plot_corner_numbers!
export colorbar!, add_colorbar
export plot_element_field!, plot_face_field!, plot_corner_field!, plot_averaged_velocity_field!
export element_field_from_function, face_field_from_functions, corner_field_from_function
export L1_norm_face, L1_norm_element, L2_norm_corner, L2_norm_face, L2_norm_element
export face_field_average_to_element_field
export get_corner_coordinates, get_face_coordinates, get_element_coordinates, get_element_corner_coordinates
export get_element_size, get_face_size
export transform_p4est_coordinates, transform_p4est_level_to_element_size
export hanging_corner_get_big_element, c2f_small, big_face_get_small_faces, face_get_big_face

export CORNERS_PER_ELEMENT, ELEMENTS_PER_CORNER, FACES_PER_ELEMENT
export LEFT, RIGHT, UP, DOWN, DOWN_LEFT, DOWN_RIGHT, UP_LEFT, UP_RIGHT
export CORNER_HANG_LEFT, CORNER_HANG_RIGHT, CORNER_HANG_DOWN, CORNER_HANG_UP
export CORNER_BOUNDARY_LEFT_MASK, CORNER_BOUNDARY_RIGHT_MASK, CORNER_BOUNDARY_DOWN_MASK, CORNER_BOUNDARY_UP_MASK
export CORNERS_PER_FACE, FACES_PER_CORNER

export CornerType, corner_regular, corner_hanging, corner_boundary
export FaceType, face_regular, face_boundary, face_big, face_first, face_second
export FaceDirection, face_vx, face_vy

export assemble_system, sol2vp, vp2sol
export compute_Txy
export StokesBoundaryConditions, stokes_boundary_free_slip, stokes_boundary_dirichlet
export populate_stokes_boundary_velocities

export  locate_point, get_point_velocity

# Note: there are many usages of @assert. For maximum performance, remove them.

const CORNERS_PER_ELEMENT=4
const ELEMENTS_PER_CORNER=4
const DOWN_LEFT=1
const DOWN_RIGHT=2
const UP_LEFT=3
const UP_RIGHT=4

function opp_diag(d)
  if d == DOWN_LEFT; return UP_RIGHT; end
  if d == DOWN_RIGHT; return UP_LEFT; end
  if d == UP_LEFT; return DOWN_RIGHT; end
  if d == UP_RIGHT; return DOWN_LEFT; end
  @assert false
end

const FACES_PER_ELEMENT=4
const ELEMENTS_PER_FACE=2
const LEFT=1
const RIGHT=2
const DOWN=3
const UP=4

const CORNERS_PER_FACE=2
const FACES_PER_CORNER=4

struct ElementInfo # the same as in a "quad" in p4est
 x::Int32  # values as in p4est (large integers)
 y::Int32  # Note: it might be better to name these e.g. "x_raw" to avoid confusion
 level::Int8
end

@enum FaceType face_regular face_boundary face_big face_first face_second
@enum FaceDirection face_vx face_vy
struct FaceInfo
 face_type::FaceType
 direction::FaceDirection  # FIXME: naming inconsistent - should be face_direction, probably
end

@enum CornerType corner_undefined corner_regular corner_hanging corner_boundary
struct CornerInfo
 corner_type::CornerType
 info::Int8        # Each bit is either 0 or 1 for the relative level (usual Z ordering)
                   # The value 0b1111 is not used in the regular case (use 0b0000)
                   # for hanging nodes, there are only 4 possible values
                   # for boundary nodes,
                   #   if the first two bits are both 1, this means corner, and
                   #     the second two bits tell you which of the 4 corners (z order: dl dr ul ur)
                   #   otherwise, the first two bits give you the levels of the two neighbors and
                   #      the second two bits tell you which of 4 faces (left right down up)
end

# Convenient constants for the info field
const CORNER_HANG_LEFT=0b0101
const CORNER_HANG_RIGHT=0b1010
const CORNER_HANG_UP=0b1100
const CORNER_HANG_DOWN=0b0011

const CORNER_BOUNDARY_DOWN_LEFT  = 0b1100
const CORNER_BOUNDARY_DOWN_RIGHT = 0b1101
const CORNER_BOUNDARY_UP_LEFT    = 0b1110
const CORNER_BOUNDARY_UP_RIGHT   = 0b1111

# use these and extract the levels from the first two bits
const CORNER_BOUNDARY_LEFT_MASK  = 0b0000
const CORNER_BOUNDARY_RIGHT_MASK = 0b0001
const CORNER_BOUNDARY_DOWN_MASK  = 0b0010
const CORNER_BOUNDARY_UP_MASK    = 0b0011

mutable struct SimpleTreeMesh
  nc::Int # The number of corners, including hanging corners
  nf::Int # The number of faces, with *three* entries for each split face.
          # Ordering is always big, first small, second small
  ne::Int # The number of elements
  e::Vector{ElementInfo}
  c::Vector{CornerInfo}
  f::Vector{FaceInfo}
  e2c # column i gives the corners which touch each element.
  c2e
  e2f
  f2e
  f2c
  c2f
  coordinate_scale::Float64
  coordinate_offset_x::Float64
  coordinate_offset_y::Float64
  SimpleTreeMesh() = new()
end

include("SimpleTreeMeshes_creation.jl")
include("SimpleTreeMeshes_functions.jl")
include("SimpleTreeMeshes_particles.jl")
include("SimpleTreeMeshes_plot.jl")
include("SimpleTreeMeshes_field.jl")
include("SimpleTreeMeshes_stokes.jl")

end  # module
