export set_coords
function set_coords(c::GraphicsContext, x, y, w, h, l, r, t, b)
    Base.depwarn("set_coords is deprecated; use set_coordinates or inner_canvas instead", :set_coords)
    inner_canvas(c::GraphicsContext, x, y, w, h, l, r, t, b)
end
function set_coords(c::GraphicsContext, device::BoundingBox, user::BoundingBox)
    Base.depwarn("set_coords is deprecated; use set_coordinates or inner_canvas instead", :set_coords)
    inner_canvas(c::GraphicsContext, device, user)
end
function set_coords(c::GraphicsContext, user::BoundingBox)
    Base.depwarn("set_coords is deprecated; use set_coordinates instead", :set_coords)
    inner_canvas(c::GraphicsContext, BoundingBox(0, width(c), 0, height(c)), user)
end
