# DCP attribute generic methods
setGeneric("is_zero", function(object) { standardGeneric("is_zero") })
setGeneric("is_positive", function(object) { standardGeneric("is_positive") })
setGeneric("is_negative", function(object) { standardGeneric("is_negative") })
setGeneric("is_unknown", function(object) { standardGeneric("is_unknown") })

setGeneric("is_constant", function(object) { standardGeneric("is_constant") })
setGeneric("is_affine", function(object) { standardGeneric("is_affine") })
setGeneric("is_convex", function(object) { standardGeneric("is_convex") })
setGeneric("is_concave", function(object) { standardGeneric("is_concave") })
setGeneric("is_dcp", function(object) { standardGeneric("is_dcp") })

setGeneric("size", function(object) { standardGeneric("size") })
setGeneric("dcp_curvature", function(monotonicity, func_curvature, arg_sign, arg_curvature) { standardGeneric("dcp_curvature") })
setGeneric("DCPAttr.mul_elemwise", function(lh_exp, rh_exp) { standardGeneric("DCPAttr.mul_elemwise") })

# Expression generic methods
setGeneric("value", function(object) { standardGeneric("value") })
setGeneric("get_data", function(object) { standardGeneric("get_data") })
setGeneric("init_dcp_attr", function(object) { standardGeneric("init_dcp_attr") })
setGeneric("curvature", function(object) { standardGeneric("curvature") })
setGeneric("is_scalar", function(object) { standardGeneric("is_scalar") })
setGeneric("is_vector", function(object) { standardGeneric("is_vector") })
setGeneric("is_matrix", function(object) { standardGeneric("is_matrix") })

setGeneric("variables", function(object) { standardGeneric("variables") })
setGeneric("parameters", function(object) { standardGeneric("parameters") })
setGeneric("canonical_form", function(object) { standardGeneric("canonical_form") })
setGeneric("canonicalize", function(object) { standardGeneric("canonicalize") })

# Atom generic methods
setGeneric("validate_args", function(object) { standardGeneric("validate_args") })
setGeneric("shape_from_args", function(object) { standardGeneric("shape_from_args") })
setGeneric("sign_from_args", function(object) { standardGeneric("sign_from_args") })
setGeneric("func_curvature", function(object) { standardGeneric("func_curvature") })
setGeneric("monotonicity", function(object) { standardGeneric("monotonicity") })
setGeneric("get_data", function(object) { standardGeneric("get_data") })
setGeneric("name", function(object) { standardGeneric("name") })

setGeneric("Atom.dcp_curvature", function(curvature, args, monotonicities) { standardGeneric("Atom.dcp_curvature") })
setGeneric("graph_implementation", function(object, arg_objs, size, data) { standardGeneric("graph_implementation") })
setGeneric("sum_squares", function(expr) { standardGeneric("sum_squares") })

# Constraint generic methods
setGeneric("id", function(object) { standardGeneric("id") })
setGeneric("violation", function(object) { standardGeneric("violation") })
setGeneric("num_cones", function(object) { standardGeneric("num_cones") })
setGeneric("cone_size", function(object) { standardGeneric("cone_size") })
