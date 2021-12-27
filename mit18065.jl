### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ 23e4c789-2cc8-4565-9042-02168b16beac
begin
	using LinearAlgebra
	using PlutoUI

	function pp(xs...)
		for x ∈ xs
			if isa(x, AbstractArray)
				println(repr("text/plain", x))
			else
				print(x)
			end
		end
		println()
	end
end

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
	@assert z == zero(z)
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
md"""
``m`` is 3.

``n`` is at least 3.

rank ``r`` is 3.
"""

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
PlutoUI.with_terminal() do
	C = rand(-9:9, 3, 2)
	R = rand(-9:9, 2, 3)
	A = C * R
	A2 = kron([0 1
		       0 1], A)
    C2 = kron([1
	           1], C)
	R2 = kron([0 1], R)
	@assert A2 == C2 * R2

	pp("C2:\n", C2)
	pp("R2:\n", R2)
	pp("A2:\n", A2)
end

# ╔═╡ df291cf8-beb9-4481-bf71-0cc67a6dab7a
md"""
### 2. Multiplying and Factoring Matrices
"""

# ╔═╡ d9ca6fbe-dffb-4c5e-b1cf-0c80c7552088
md"""
**(2.2)** Suppose ``\mathbf{a}`` and ``\mathbf{b}`` are column vectors with components ``a_1,...,a_m`` and ``b_1,...,b_p``. Can you multiply ``\mathbf{a}`` times ``\mathbf{b}^\mathsf{T}`` (yes or no)? What is the shape of the answer ``\mathbf{ab}^\mathsf{T}``? What number is in row ``i``, column ``j`` of ``\mathbf{ab}^\mathsf{T}``? What can you say about ``\mathbf{aa}^\mathsf{T}``?
"""

# ╔═╡ d4f90e59-3b6b-4b51-b1ae-5c0f2b2aa8ea
PlutoUI.with_terminal() do
	m = 3
	p = 5
	a = rand(-9:9, m, 1)
	b = rand(-9:9, p, 1)
	abT = a * transpose(b)
	pp("a:\n", a)
	pp("b:\n", b)
	pp("abT:\n", abT)

	@assert size(abT) == (m, p)

	i = rand(1:m)
	j = rand(1:p)
	@assert abT[i,j] == a[i] * b[j]

	@assert issymmetric(a * transpose(a))
end

# ╔═╡ ced55756-9423-4f66-a755-fa7afb61c159
md"""
**(2.6)** If ``A`` has columns ``\mathbf{a}_1``, ``\mathbf{a}_2``, ``\mathbf{a}_3`` and ``B = I`` is the identity matrix, what are the rank one matrices ``\mathbf{a}_1\mathbf{b}^*_1`` and ``\mathbf{a}_2\mathbf{b}^*_2`` and ``\mathbf{a}_3\mathbf{b}^*_3``? They should add to ``AI = A``.
"""

# ╔═╡ 3d07a918-0263-4d1e-9187-b0adbbe20203
PlutoUI.with_terminal() do
	a1 = rand(-9:9, 5, 1)
	a2 = rand(-9:9, 5, 1)
	a3 = rand(-9:9, 5, 1)
	A = [a1 a2 a3]
	pp("A:\n", A)

	B = I(3)
	pp("B:\n", B)

	a1b1 = a1 * B[:,1]'
	a2b2 = a2 * B[:,2]'
	a3b3 = a3 * B[:,3]'
	@assert a1b1 + a2b2 + a3b3 == A

	pp("a1:\n", a1)
	pp("b'1:\n", B[:,1]')
	pp("a1b'1:\n", a1b1)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.27"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "d7fa6237da8004be601e19bd6666083056649918"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "fed057115644d04fba7f4d768faeeeff6ad11a60"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.27"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
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
# ╟─5af2ab12-84bb-4cf5-b632-ffdfce21b8d0
# ╟─a3511d6e-c7dc-44cc-a6b3-43e7440e32cc
# ╠═005cde91-f411-40e8-a903-01506402c2be
# ╟─df291cf8-beb9-4481-bf71-0cc67a6dab7a
# ╟─d9ca6fbe-dffb-4c5e-b1cf-0c80c7552088
# ╠═d4f90e59-3b6b-4b51-b1ae-5c0f2b2aa8ea
# ╟─ced55756-9423-4f66-a755-fa7afb61c159
# ╠═3d07a918-0263-4d1e-9187-b0adbbe20203
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
