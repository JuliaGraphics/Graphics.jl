var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference-1","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/#Geometric-primitives-1","page":"Reference","title":"Geometric primitives","text":"","category":"section"},{"location":"reference/#","page":"Reference","title":"Reference","text":"Vec2\nPoint\nBoundingBox","category":"page"},{"location":"reference/#Graphics.Vec2","page":"Reference","title":"Graphics.Vec2","text":"Vec2(x, y) -> v\n\nCreate a Cartesian representation v of a vector (or point) in two dimensions.\n\n\n\n\n\n","category":"type"},{"location":"reference/#Graphics.Point","page":"Reference","title":"Graphics.Point","text":"Point(x, y) -> p\n\nCreate a Cartesian representation p of a point in two dimensions. Point is an alias of Vec2.\n\n\n\n\n\n","category":"type"},{"location":"reference/#Graphics.BoundingBox","page":"Reference","title":"Graphics.BoundingBox","text":"BoundingBox(xmin, xmax, ymin, ymax) -> bb\n\nCreate a representation bb of a rectangular region, specifying the coordinates of the horizontal (x) and vertical (y) edges.\n\n\n\n\n\n","category":"type"},{"location":"reference/#Geometry-API-1","page":"Reference","title":"Geometry API","text":"","category":"section"},{"location":"reference/#","page":"Reference","title":"Reference","text":"aspect_ratio\ncenter\ndeform\ndiagonal\nisinside\nshift\nheight\nwidth\nxmin\nxmax\nymin\nymax\nxrange\nyrange","category":"page"},{"location":"reference/#Graphics.aspect_ratio","page":"Reference","title":"Graphics.aspect_ratio","text":"aspect_ratio(bb::BoundingBox) -> r\n\nCompute the ratio r of the height and width of bb.\n\nExample\n\njulia> bb = BoundingBox(Point(0, 0), Point(1920, 1080)); # landscape\n\njulia> aspect_ratio(bb)\n0.5625\n\njulia> rationalize(ans)\n9//16\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.center","page":"Reference","title":"Graphics.center","text":"center(obj) -> p::Point\n\nCompute the center coordinate of obj.\n\nnote: Note\nThe fallback implementation of this function returns the center of bounding box, not the geometric center, or centroid.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.deform","page":"Reference","title":"Graphics.deform","text":"deform(bb::BoundingBox, Δl, Δr, Δt, Δb) -> bbnew\n\nAdd Δl (left), Δr (right), Δt (top), and Δb (bottom) to the edges of a BoundingBox. The sign of each value follows the positive direction of the axis. The \"top\" and \"bottom\" are representations when the y-axis is directed downward.\n\nExample\n\njulia> bb = BoundingBox(Point(1, 1), Point(10, 10));\n\njulia> bbnew = deform(bb, 0.5, -1.0, -0.25, 1.0);\n\njulia> xrange(bbnew), yrange(bbnew)\n((1.5, 9.0), (0.75, 11.0))\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.diagonal","page":"Reference","title":"Graphics.diagonal","text":"diagonal(obj) -> diag\n\nGet the diagonal length of obj. The fallback implementation of this function returns the diagonal length of the bounding box of obj.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.isinside","page":"Reference","title":"Graphics.isinside","text":"isinside(bb::BoundingBox, p::Point) -> tf::Bool\nisinside(bb::BoundingBox, x, y) -> tf::Bool\n\nDetermine whether the point lies within bb.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.shift","page":"Reference","title":"Graphics.shift","text":"shift(bb::BoundingBox, Δx, Δy) -> bbnew\n\nShift center by (Δx, Δy), keeping width & height fixed.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.height","page":"Reference","title":"Graphics.height","text":"hieght(obj) -> h\n\nGet the vertical length of obj.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.width","page":"Reference","title":"Graphics.width","text":"width(obj) -> w\n\nGet the horizontal length of obj.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.xmin","page":"Reference","title":"Graphics.xmin","text":"xmin(obj) -> xmin\n\nGet the minimum x coordinate of the bounding box of obj.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.xmax","page":"Reference","title":"Graphics.xmax","text":"xmax(obj) -> xmax\n\nGet the maximum x coordinate of the bounding box of obj.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.ymin","page":"Reference","title":"Graphics.ymin","text":"ymin(obj) -> ymin\n\nGet the minimum x coordinate of the bounding box of obj.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.ymax","page":"Reference","title":"Graphics.ymax","text":"ymax(obj) -> ymax\n\nGet the maximum y coordinate of the bounding box of obj.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.xrange","page":"Reference","title":"Graphics.xrange","text":"xrange(obj) -> (xmin, xmax)\n\nGet the horizontal range of the bounding box that minimally contains obj.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.yrange","page":"Reference","title":"Graphics.yrange","text":"yrange(obj) -> (ymin, ymax)\n\nGet the vertical range of the bounding box that minimally contains obj.\n\n\n\n\n\n","category":"function"},{"location":"reference/#d-drawing-contexts-1","page":"Reference","title":"2d drawing contexts","text":"","category":"section"},{"location":"reference/#","page":"Reference","title":"Reference","text":"GraphicsDevice\nGraphicsContext\ncreategc\ngetgc","category":"page"},{"location":"reference/#Graphics.GraphicsDevice","page":"Reference","title":"Graphics.GraphicsDevice","text":"GraphicDevice\n\nAn abstract graphics output device; can create GraphicsContexts.\n\n\n\n\n\n","category":"type"},{"location":"reference/#Graphics.GraphicsContext","page":"Reference","title":"Graphics.GraphicsContext","text":"GraphicsContext\n\nAn abstract object that can actually be drawn to.\n\n\n\n\n\n","category":"type"},{"location":"reference/#Graphics.creategc","page":"Reference","title":"Graphics.creategc","text":"creategc(gd::GraphicsDevice) -> gc::GraphicContext\n\nCreate a new GraphicContext.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.getgc","page":"Reference","title":"Graphics.getgc","text":"getgc(obj) -> gc::GraphicContext\n\nGet a GraphicsContext from something that might be drawable.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Coordinate-systems-1","page":"Reference","title":"Coordinate systems","text":"","category":"section"},{"location":"reference/#","page":"Reference","title":"Reference","text":"set_coordinates\nreset_transform\nrotate\nscale\ntranslate\nuser_to_device!\ndevice_to_user!\nuser_to_device_distance!\ndevice_to_user_distance!\nuser_to_device\ndevice_to_user","category":"page"},{"location":"reference/#Graphics.set_coordinates","page":"Reference","title":"Graphics.set_coordinates","text":"set_coordinates(gc::GraphicsContext, device::BoundingBox, user::BoundingBox)\nset_coordinates(gc::GraphicsContext, user::BoundingBox)\n\nSet the device->user coordinate transformation of c so that device, expressed in \"device coordinates\" (pixels), is equivalent to user as expressed in \"user coordinates\". If device is omitted, it defaults to the full span of gc, BoundingBox(0, width(gc), 0, height(gc)).\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.reset_transform","page":"Reference","title":"Graphics.reset_transform","text":"reset_transform(gc::GraphicsContext)\n\nReset the current transformation.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.rotate","page":"Reference","title":"Graphics.rotate","text":"rotate(p::Vec2, angle::Real, o::Vec2 = Vec2(0, 0)) -> pnew::Vec2\n\nRotate p around o by angle (in radians).\n\nnote: Note\nThe direction of rotation for positive angles is from the positive x-axis toward the positive y-axis. For example, the direction of rotation is \"clockwise\" when the x-axis is directed right and the y-axis is directed downward.\n\nExample\n\njulia> rotate(Vec2(2, 1), 0.5π, Vec2(1, 1))\nVec2(1.0, 2.0)\n\n\n\n\n\nrotate(bb::BoundingBox, angle::Real, o::Point) -> bbnew\n\nRotate bb around o by angle (in radians), returning the BoundingBox that encloses the vertices of the rotated box.\n\n\n\n\n\nrotate(gc::GraphicsContext, angle)\n\nRotate the user-space axes by angle (in radians). The rotation takes places after any existing transformation.\n\nSee also: rotate(p::Vec2, angle::Real, o::Vec2).\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.scale","page":"Reference","title":"Graphics.scale","text":"rotate(gc::GraphicsContext, sx, sy)\n\nScale the user-space x-axis and y-axis by sx and sy respectively. The scaling takes places after any existing transformation.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.translate","page":"Reference","title":"Graphics.translate","text":"translate(gc::GraphicsContext, Δx, Δy)\n\nTranslate the user-space origin by (Δx, Δy). The translation takes places after any existing transformation.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.user_to_device!","page":"Reference","title":"Graphics.user_to_device!","text":"user_to_device!(gc::GraphicsContext, c::Vector{Float64})\n\nTransform a coordinate c from the user space to the device space.\n\nSee also: user_to_device, device_to_user!\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.device_to_user!","page":"Reference","title":"Graphics.device_to_user!","text":"device_to_user!(gc::GraphicsContext, c::Vector{Float64})\n\nTransform a coordinate c from the device space to the user space.\n\nSee also: device_to_user, user_to_device!\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.user_to_device_distance!","page":"Reference","title":"Graphics.user_to_device_distance!","text":"user_to_device_distance!(gc::GraphicsContext, d::Vector{Float64})\n\nTransform a distance vector d from the user space to the device space. This function is similar to the device_to_user! except that the translation components will be cancelled.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.device_to_user_distance!","page":"Reference","title":"Graphics.device_to_user_distance!","text":"device_to_user_distance!(gc::GraphicsContext, d::Vector{Float64})\n\nTransform a distance vector d from the device space to the user space. This function is similar to the user_to_device! except that the translation components will be cancelled.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.user_to_device","page":"Reference","title":"Graphics.user_to_device","text":"user_to_device(gc::GraphicsContext, x, y) -> (xd, yd)\n\nTransform a user space coordinate (x, y) to the device space coordinate (xd, yd).\n\nSee also: user_to_device!, device_to_user\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.device_to_user","page":"Reference","title":"Graphics.device_to_user","text":"device_to_user(gc::GraphicsContext, x, y) -> (xu, yu)\n\nTransform a device space coordinate (x, y) to the user space coordinate (xu, yu).\n\nSee also: device_to_user!, user_to_device\n\n\n\n\n\n","category":"function"},{"location":"reference/#Lines-1","page":"Reference","title":"Lines","text":"","category":"section"},{"location":"reference/#","page":"Reference","title":"Reference","text":"set_line_width\nset_dash","category":"page"},{"location":"reference/#Graphics.set_line_width","page":"Reference","title":"Graphics.set_line_width","text":"set_line_width(gc::GraphicsContext, w)\n\nSet the current line width (in device-depended units). The actual width and aspect-ratio on the screen may be affected by the transformation.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.set_dash","page":"Reference","title":"Graphics.set_dash","text":"set_dash(gc::GraphicsContext, dashes::Vector{Float64}, offset)\n\nSet the dash pattern. The dashes is a Vector of positive lengths. The odd-numbered elements represent the length of the \"on\" state, and the even-numbered elements represent the length of the \"off\" (blank) state. The offset specifies an offset at which the stroke begins. If dashes is empty, dashing is disabled.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Colors-and-painting-(drawing-attributes)-1","page":"Reference","title":"Colors and painting (drawing attributes)","text":"","category":"section"},{"location":"reference/#","page":"Reference","title":"Reference","text":"set_source\nset_source_rgb\nset_source_rgba\nsave\nrestore","category":"page"},{"location":"reference/#Graphics.set_source","page":"Reference","title":"Graphics.set_source","text":"set_source(gc::GraphicsContext, src)\n\nSet the source pattern to src. If the src is a Color, it is used as the source color.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.set_source_rgb","page":"Reference","title":"Graphics.set_source_rgb","text":"set_source_rgb(gc::GraphicsContext, r, g, b)\n\nSet the source pattern to an opaque color RGB(r, g, b). The color components are in the range [0, 1], not [0, 255].\n\nSee also: set_source_rgba.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.set_source_rgba","page":"Reference","title":"Graphics.set_source_rgba","text":"set_source_rgba(gc::GraphicsContext, r, g, b, a)\n\nSet the source pattern to a transparent color RGBA(r, g, b, a). The color and alpha components are in the range [0, 1], not [0, 255].\n\nSee also: set_source_rgb.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.save","page":"Reference","title":"Graphics.save","text":"save(gc::GraphicsContext)\n\nSave the copy of current context gc. The context is saved onto an internal stack, for example.\n\nSee also: restore\n\nwarning: Warning\nThe function name save conflicts with the save in the FileIO package, which is likely to be used with the Graphics package.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.restore","page":"Reference","title":"Graphics.restore","text":"restore(gc::GraphicsContext)\n\nRestore the context saved by a preceding save. The state at the time of the call is overwritten with the restored context, and the restored context is removed from an internal stack, for example.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Clipping-1","page":"Reference","title":"Clipping","text":"","category":"section"},{"location":"reference/#","page":"Reference","title":"Reference","text":"clip\nclip_preserve\nreset_clip\ninner_canvas","category":"page"},{"location":"reference/#Graphics.clip","page":"Reference","title":"Graphics.clip","text":"clip(gc::GraphicsContext)\n\nSet a new clipping region based on the current path and fill within the context.\n\nSee also: reset_clip, clip_preserve.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.clip_preserve","page":"Reference","title":"Graphics.clip_preserve","text":"clip_preserve(gc::GraphicsContext)\n\nSet a new clipping region based on the current path and fill within the context. Unlike the clip function, this function preserves the path within the context.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.reset_clip","page":"Reference","title":"Graphics.reset_clip","text":"reset_clip(gc::GraphicsContext)\n\nRemove the clipping region set by the clip or clip_preserve function.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.inner_canvas","page":"Reference","title":"Graphics.inner_canvas","text":"inner_canvas(gc::GraphicsContext, device::BoundingBox, user::BoundingBox)\ninner_canvas(gc::GraphicsContext, x, y, w, h, l, r, t, b)\n\nCreate a rectangular drawing area inside device (represented in device-coordinates), giving it user-coordinates user. Any drawing that occurs outside this box is clipped.\n\nx, y, w, and h are an alternative parametrization of device, and l, r, t, b parametrize user.\n\nSee also: set_coordinates.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Paths-1","page":"Reference","title":"Paths","text":"","category":"section"},{"location":"reference/#","page":"Reference","title":"Reference","text":"move_to\nline_to\nrel_line_to\nrel_move_to\nnew_path\nnew_sub_path\nclose_path\narc","category":"page"},{"location":"reference/#Graphics.move_to","page":"Reference","title":"Graphics.move_to","text":"move_to(gc::GraphicsContext, x, y)\n\nBegin a new sub-path. The current point will be moved to (x, y).\n\nSee also: rel_move_to.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.line_to","page":"Reference","title":"Graphics.line_to","text":"line_to(gc::GraphicsContext, x, y)\n\nAdd a line to the current path from the current point to the position (x, y). The current point will be moved to (x, y).\n\nSee also: rel_line_to.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.rel_line_to","page":"Reference","title":"Graphics.rel_line_to","text":"rel_line_to(gc::GraphicsContext, Δx, Δy)\n\nAdd a line to the current path from the current point of (x, y) to the position (x + Δx, y + Δy). The current point will be moved to (x + Δx, y + Δy).\n\nSee also: line_to.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.rel_move_to","page":"Reference","title":"Graphics.rel_move_to","text":"rel_move_to(gc::GraphicsContext, Δx, Δy)\n\nBegin a new sub-path. The current point will be moved to (x + Δx, y + Δy), where (x, y) is the previous current point.\n\nSee also: move_to.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.new_path","page":"Reference","title":"Graphics.new_path","text":"new_path(gc::GraphicsContext)\n\nClear the current path.\n\nnote: Note\nDepending on the backend, the new path may actually begin when a path element is added.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.new_sub_path","page":"Reference","title":"Graphics.new_sub_path","text":"new_sub_path(gc::GraphicsContext)\n\nBegin a new sub-path. The current point will be cleared.\n\nSee also: move_to.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.close_path","page":"Reference","title":"Graphics.close_path","text":"close_path(gc::GraphicsContext)\n\nAdd a line to the current path from the current point to the beginning of the current sub-path and closes the sub-path. The current point will be changed to the joined point.\n\nThere is a difference between closing a subpath and drawing a line to the equivalent coordinate. This difference might be visualized as a difference in drawing stroke endpoints.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.arc","page":"Reference","title":"Graphics.arc","text":"arc(gc::GraphicsContext, xc, yc, radius, angle1, angle2)\n\nAdd a circular arc with the specified radius to the current path. The arc is centered at (xc, yc), begins at angle1 and ends at angle2. The angle1 and angle2 are in radians. The arc will be drawn in the direction of increasing angles.\n\n\n\n\n\n","category":"function"},{"location":"reference/#High-level-paths-1","page":"Reference","title":"High-level paths","text":"","category":"section"},{"location":"reference/#","page":"Reference","title":"Reference","text":"rectangle\ncircle\npolygon","category":"page"},{"location":"reference/#Graphics.rectangle","page":"Reference","title":"Graphics.rectangle","text":"rectangle(gc::GraphicsContext, x, y, width, height)\nrectangle(gc::GraphicsContext, user::BoundingBox)\n\nAdd a sub-path rectangle to the current path. The x and y specify the (typically the upper left) corner coordinate of the rectangle, and the width and height specify the size.\n\nYou can also specify the position and size by a BoundingBox in user-space coordinate.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.circle","page":"Reference","title":"Graphics.circle","text":"circle(ctx::GraphicsContext, x, y, r)\n\nAdd a sub-path circle to the current path. The x and y specify the center coordinate of the circle, and the r specifies the radius.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.polygon","page":"Reference","title":"Graphics.polygon","text":"polygon(gc::GraphicsContext, verts::Matrix, idx::Vector)\n\nAdd a closed sub-path polygon with the given vertices. The verts is a collection of vertex coordinates in the following matrix form:\n\n[x1 x2 x3 ... xn;\n y1 y2 y3 ... yn]\n\nThe idx is a vector of vertex indices, i.e. the matrix column numbers.\n\ntip: Tip\nYou can reuse the vertex coordinates by specifying the same index in idx. This is useful when drawing meshes.\n\n\n\n\n\npolygon(gc::GraphicsContext, points::AbstractVector)\n\nAdd a closed sub-path polygon with the given vertices to .\n\n\n\n\n\n","category":"function"},{"location":"reference/#Fill-and-stroke-1","page":"Reference","title":"Fill and stroke","text":"","category":"section"},{"location":"reference/#","page":"Reference","title":"Reference","text":"fill\nfill_preserve\npaint\nstroke\nstroke_preserve\nstroke_transformed\nstroke_transformed_preserve","category":"page"},{"location":"reference/#Base.fill","page":"Reference","title":"Base.fill","text":"fill(gc::GraphicsContext)\n\nFill the current path according to the current fill rule. The current path will be cleared from the context.\n\nSee also: fill_preserve.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.fill_preserve","page":"Reference","title":"Graphics.fill_preserve","text":"fill_preserve(gc::GraphicsContext)\n\nFill the current path according to the current fill rule. Unlike the fill function, this function preserves the current path within the context.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.paint","page":"Reference","title":"Graphics.paint","text":"paint(gc::GraphicsContext)\n\nPaint the current source everywhere within the current clipping region.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.stroke","page":"Reference","title":"Graphics.stroke","text":"stroke(gc::GraphicsContext)\n\nStroke the current path according to the current stroke style. The current path will be cleared from the context.\n\nnote: Note\nThe stroke function ignores the current transformation. If you want to apply the current transformation to the stroke, use stroke_transformed.\n\nSee also: stroke_preserve.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.stroke_preserve","page":"Reference","title":"Graphics.stroke_preserve","text":"stroke_preserve(gc::GraphicsContext)\n\nStroke the current path according to the current stroke style. Unlike the stroke function, this function preserves the current path within the context.\n\nnote: Note\nThe stroke_preserve function ignores the current transformation. If you want to apply the current transformation to the stroke, use stroke_transformed_preserve.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.stroke_transformed","page":"Reference","title":"Graphics.stroke_transformed","text":"stroke_transformed(gc::GraphicsContext)\n\nStroke the current path according to the current stroke style and the current transformation.\n\nSee also: stroke.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Graphics.stroke_transformed_preserve","page":"Reference","title":"Graphics.stroke_transformed_preserve","text":"stroke_transformed_preserve(gc::GraphicsContext)\n\nStroke the current path according to the current stroke style and the current transformation. Unlike the stroke_transformed function, this function preserves the current path within the context.\n\nSee also: stroke_preserve.\n\n\n\n\n\n","category":"function"},{"location":"#Graphics.jl-1","page":"Graphics.jl","title":"Graphics.jl","text":"","category":"section"},{"location":"#","page":"Graphics.jl","title":"Graphics.jl","text":"Graphics.jl is an abstraction layer for graphical operations in Julia. Its goal is to allow developers to write graphical programs in a manner independent of the particular graphical backend. One needs to load a backend package that implements the operations in its API; currently, Cairo.jl is the only such backend.","category":"page"},{"location":"#","page":"Graphics.jl","title":"Graphics.jl","text":"To get an organized overview of the API, try typing ?Graphics at the Julia REPL. You can see the same information in greater detail on the next page.","category":"page"}]
}
