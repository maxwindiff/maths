### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ 23e4c789-2cc8-4565-9042-02168b16beac
using LinearAlgebra

# ╔═╡ 246ffefa-6614-11ec-09c8-1b72223a1641
html"""
	<div>Happy holiday! Remember to take care of yourself and your loved ones!</div>
<div id="snow"></div>
<style>
	body:not(.disable_ui):not(.more-specificity) {
        background-color:#e9ecff;
    }
	pluto-output{
		border-radius: 0px 8px 0px 0px;
        background-color:#e9ecff;
	}
	#snow {
        position: fixed;
    	top: 0;
    	left: 0;
    	right: 0;
    	bottom: 0;
    	pointer-events: none;
    	z-index: 1000;
	}
</style>
<script src="https://cdn.jsdelivr.net/particles.js/2.0.0/particles.min.js"></script>
<script>
        setTimeout(() => window.particlesJS("snow", {
            "particles": {
                "number": {
                    "value": 70,
                    "density": {
                        "enable": true,
                        "value_area": 800
                    }
                },
                "color": {
                    "value": "#ffffff"
                },
                "opacity": {
                    "value": 0.7,
                    "random": false,
                    "anim": {
                        "enable": false
                    }
                },
                "size": {
                    "value": 5,
                    "random": true,
                    "anim": {
                        "enable": false
                    }
                },
                "line_linked": {
                    "enable": false
                },
                "move": {
                    "enable": true,
                    "speed": 5,
                    "direction": "bottom",
                    "random": true,
                    "straight": false,
                    "out_mode": "out",
                    "bounce": false,
                    "attract": {
                        "enable": true,
                        "rotateX": 300,
                        "rotateY": 1200
                    }
                }
            },
            "interactivity": {
                "events": {
                    "onhover": {
                        "enable": false
                    },
                    "onclick": {
                        "enable": false
                    },
                    "resize": false
                }
            },
            "retina_detect": true
        }), 3000);
	</script>
"""


# ╔═╡ 3d9413a9-32aa-439e-8a5d-aa60ab15a651
md"""
### 1. The Column Space of ``A`` Contains All Vectors ``Ax``
"""

# ╔═╡ f5334940-2e47-4301-9508-ecfc7108f334
md"""
**(1.1)** Give an example where a combination of three nonzero vectors in ``\mathbf{R}^4`` is the zero vector. Then write your example in the form ``A\mathbf{x} = \mathbf{0}``. What are the shapes of ``A`` and ``\mathbf{x}`` and ``\mathbf{0}``?
"""

# ╔═╡ 246fff04-6614-11ec-22e1-d995b65a10a2
let
	A = [1 0 1
		 0 1 1
		 0 0 0
		 0 0 0]
	x = [1
	     1
	    -1]
	z = A * x
	@assert size(A) == (4,3)
	@assert size(x) == (3,)
	@assert size(z) == (4,)
	z
end

# ╔═╡ 96c73c4a-f62c-41d6-bcc0-13944d3fd3c7
md"""
**(1.4)** Suppose ``A`` is the 3 by 3 matrix ``\mathbf{ones}(3, 3)`` of all ones. Find two independent vectors ``\mathbf{x}`` and ``\mathbf{y}`` that solve ``A\mathbf{x}=\mathbf{0}`` and ``A\mathbf{y}=\mathbf{0}``. Write that first equation ``A\mathbf{x}=\mathbf{0}`` (with numbers) as a combination of the columns of ``A``. Why don’t I ask for a third independent vector with ``A\mathbf{z} = \mathbf{0}``?
"""

# ╔═╡ 6f2f5e4e-7843-46b5-8647-3f681cfbbae0
let
	A = ones(3, 3)
	x = [1
		-1
		0]
	@assert A * x == zero(A * x)
	y = [1
		 0
	    -1]
	@assert A * y == zero(A * y)
	# A has rank 1, so dimension of null space is 3 - 1 = 2
	@assert rank(A) == 1
	@assert size(nullspace(A), 2) == 2
end

# ╔═╡ b344b879-d30a-420d-9e8c-410ab75c75a5
md"""
**(1.9)** Suppose the column space of an ``m`` by ``n`` matrix is all of ``\mathbf{R}^3``. What can you say about ``m``? What can you say about ``n``? What can you say about the rank ``r``?
"""

# ╔═╡ 5af2ab12-84bb-4cf5-b632-ffdfce21b8d0
let
	# m is at 3, n is at least 3, rank is 3
end

# ╔═╡ a3511d6e-c7dc-44cc-a6b3-43e7440e32cc
md"""
**(1.18)** If ``A = CR``, what are the ``CR`` factors of the matrix
```math
\begin{equation*}
\begin{bmatrix}
0 & A \\
0 & A \\
\end{bmatrix}
\end{equation*}
```
"""

# ╔═╡ 005cde91-f411-40e8-a903-01506402c2be
let
	C = rand(-5:5, 3, 2)
	R = rand(-5:5, 2, 3)
	A = C * R
	A2 = kron([0 1
		       0 1], A)
    C2 = kron([1
	           1], C)
	R2 = kron([0 1], R)
	@assert A2 == C2 * R2
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╟─246ffefa-6614-11ec-09c8-1b72223a1641
# ╠═23e4c789-2cc8-4565-9042-02168b16beac
# ╟─3d9413a9-32aa-439e-8a5d-aa60ab15a651
# ╟─f5334940-2e47-4301-9508-ecfc7108f334
# ╠═246fff04-6614-11ec-22e1-d995b65a10a2
# ╟─96c73c4a-f62c-41d6-bcc0-13944d3fd3c7
# ╠═6f2f5e4e-7843-46b5-8647-3f681cfbbae0
# ╟─b344b879-d30a-420d-9e8c-410ab75c75a5
# ╠═5af2ab12-84bb-4cf5-b632-ffdfce21b8d0
# ╟─a3511d6e-c7dc-44cc-a6b3-43e7440e32cc
# ╠═005cde91-f411-40e8-a903-01506402c2be
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
