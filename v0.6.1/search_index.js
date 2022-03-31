var documenterSearchIndex = {"docs":
[{"location":"VonMises/#von-Mises-model","page":"von Mises model","title":"von Mises model","text":"","category":"section"},{"location":"VonMises/","page":"von Mises model","title":"von Mises model","text":"VonMises","category":"page"},{"location":"VonMises/#MaterialModels.VonMises","page":"von Mises model","title":"MaterialModels.VonMises","text":"VonMises(::ElasticModel; q_y)\n\nYield function\n\nf = q - q_mathrmy\n\nPlastic flow\n\ng = q - q_mathrmy\n\n\n\n\n\nVonMises(::ElasticModel, mohr_coulomb_type; c)\n\nParameters\n\nmohr_coulomb_type: :planestrain\nc: cohesion\n\n\n\n\n\n","category":"type"},{"location":"VonMises/#Methods","page":"von Mises model","title":"Methods","text":"","category":"section"},{"location":"VonMises/","page":"von Mises model","title":"von Mises model","text":"Modules = [MaterialModels]\nOrder = [:function]\nPages = [\"VonMises.jl\"]","category":"page"},{"location":"VonMises/#MaterialModels.matcalc__plasticflow__σ-Tuple{VonMises, SymmetricSecondOrderTensor{3}}","page":"von Mises model","title":"MaterialModels.matcalc__plasticflow__σ","text":"@matcalc(:plasticflow, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})\n\nCompute the plastic flow (the gradient of plastic potential function in stress space, i.e., partialgpartialsigma).\n\n\n\n\n\n","category":"method"},{"location":"VonMises/#MaterialModels.matcalc__stress__dϵ__σ-Tuple{VonMises, SymmetricSecondOrderTensor{3}, SymmetricSecondOrderTensor{3}}","page":"von Mises model","title":"MaterialModels.matcalc__stress__dϵ__σ","text":"@matcalc(:stress, model::VonMises{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})\n\nCompute the stress.\n\n\n\n\n\n","category":"method"},{"location":"VonMises/#MaterialModels.matcalc__stressall__dϵ__σ-Tuple{VonMises{LinearElastic}, SymmetricSecondOrderTensor{3}, SymmetricSecondOrderTensor{3}}","page":"von Mises model","title":"MaterialModels.matcalc__stressall__dϵ__σ","text":"@matcalc(:stressall, model::VonMises{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})\n\nCompute the stress and related variables as NamedTuple.\n\n\n\n\n\n","category":"method"},{"location":"VonMises/#MaterialModels.matcalc__yieldfunction__σ-Tuple{VonMises, SymmetricSecondOrderTensor{3}}","page":"von Mises model","title":"MaterialModels.matcalc__yieldfunction__σ","text":"@matcalc(:yieldfunction, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})\n\nCompute the yield function.\n\n\n\n\n\n","category":"method"},{"location":"WaterEOS/#Equations-of-state-of-water","page":"Equations of state of water","title":"Equations of state of water","text":"","category":"section"},{"location":"WaterEOS/#Monaghan","page":"Equations of state of water","title":"Monaghan","text":"","category":"section"},{"location":"WaterEOS/","page":"Equations of state of water","title":"Equations of state of water","text":"MonaghanWaterEOS","category":"page"},{"location":"WaterEOS/#MaterialModels.MonaghanWaterEOS","page":"Equations of state of water","title":"MaterialModels.MonaghanWaterEOS","text":"MonaghanWaterEOS(; ρ_ref, B, γ = 7)\n\nEquation\n\np = B left( left( fracrhorho_mathrmref right)^gamma -1 right)\n\nParameters\n\nρ_ref: reference density\n\n\n\n\n\n","category":"type"},{"location":"WaterEOS/#Methods","page":"Equations of state of water","title":"Methods","text":"","category":"section"},{"location":"WaterEOS/","page":"Equations of state of water","title":"Equations of state of water","text":"MaterialModels.matcalc__pressure__ρ(::MonaghanWaterEOS, ::Real)\nMaterialModels.matcalc__density__p(::MonaghanWaterEOS, ::Real)","category":"page"},{"location":"WaterEOS/#MaterialModels.matcalc__pressure__ρ-Tuple{MonaghanWaterEOS, Real}","page":"Equations of state of water","title":"MaterialModels.matcalc__pressure__ρ","text":"@matcalc(:pressure, model::MonaghanWaterEOS; ρ::Real)\n\nCompute the pressure.\n\n\n\n\n\n","category":"method"},{"location":"WaterEOS/#MaterialModels.matcalc__density__p-Tuple{MonaghanWaterEOS, Real}","page":"Equations of state of water","title":"MaterialModels.matcalc__density__p","text":"@matcalc(:density, model::MonaghanWaterEOS; p::Real)\n\nCompute the mass density.\n\n\n\n\n\n","category":"method"},{"location":"WaterEOS/#Morris","page":"Equations of state of water","title":"Morris","text":"","category":"section"},{"location":"WaterEOS/","page":"Equations of state of water","title":"Equations of state of water","text":"MorrisWaterEOS","category":"page"},{"location":"WaterEOS/#MaterialModels.MorrisWaterEOS","page":"Equations of state of water","title":"MaterialModels.MorrisWaterEOS","text":"MorrisWaterEOS(; ρ_ref, c)\n\nEquation\n\np = c^2 (rho - rho_mathrmref)\n\nParameters\n\nρ_ref: reference density\nc: speed of sound\n\n\n\n\n\n","category":"type"},{"location":"WaterEOS/#Methods-2","page":"Equations of state of water","title":"Methods","text":"","category":"section"},{"location":"WaterEOS/","page":"Equations of state of water","title":"Equations of state of water","text":"MaterialModels.matcalc__pressure__ρ(::MorrisWaterEOS, ::Real)\nMaterialModels.matcalc__density__p(::MorrisWaterEOS, ::Real)","category":"page"},{"location":"WaterEOS/#MaterialModels.matcalc__pressure__ρ-Tuple{MorrisWaterEOS, Real}","page":"Equations of state of water","title":"MaterialModels.matcalc__pressure__ρ","text":"@matcalc(:pressure, model::MorrisWaterEOS; ρ::Real)\n\nCompute the pressure.\n\n\n\n\n\n","category":"method"},{"location":"WaterEOS/#MaterialModels.matcalc__density__p-Tuple{MorrisWaterEOS, Real}","page":"Equations of state of water","title":"MaterialModels.matcalc__density__p","text":"@matcalc(:density, model::MorrisWaterEOS; p::Real)\n\nCompute the density.\n\n\n\n\n\n","category":"method"},{"location":"FluidModel/#Models-for-fluids","page":"Models for fluids","title":"Models for fluids","text":"","category":"section"},{"location":"FluidModel/#Newtonian-fluid","page":"Models for fluids","title":"Newtonian fluid","text":"","category":"section"},{"location":"FluidModel/","page":"Models for fluids","title":"Models for fluids","text":"NewtonianFluid","category":"page"},{"location":"FluidModel/#MaterialModels.NewtonianFluid","page":"Models for fluids","title":"MaterialModels.NewtonianFluid","text":"NewtonianFluid(::WaterEOS; μ, λ = -2μ/3)\n\nEquation\n\nsigma_ij = -p delta_ij + lambda d_kk delta_ij + 2mu d_ij\n\nwhere p is the pressure and d_ij is the rate of deformation tensor.\n\nParameters\n\nμ: dynamic viscosity\nλ: second coefficient of viscosity\n\n\n\n\n\n","category":"type"},{"location":"FluidModel/#Methods","page":"Models for fluids","title":"Methods","text":"","category":"section"},{"location":"FluidModel/","page":"Models for fluids","title":"Models for fluids","text":"MaterialModels.matcalc__pressure__ρ(::NewtonianFluid, ::Real)\nMaterialModels.matcalc__density__p(::NewtonianFluid, ::Real)\nMaterialModels.matcalc__stress__d__ρ(::NewtonianFluid, ::SymmetricSecondOrderTensor{3}, ::Real)","category":"page"},{"location":"FluidModel/#MaterialModels.matcalc__pressure__ρ-Tuple{NewtonianFluid, Real}","page":"Models for fluids","title":"MaterialModels.matcalc__pressure__ρ","text":"@matcalc(:pressure, model::NewtonianFluid{EOS}; ρ::Real)\n\nCompute the pressure based on the equation of state EOS.\n\n\n\n\n\n","category":"method"},{"location":"FluidModel/#MaterialModels.matcalc__density__p-Tuple{NewtonianFluid, Real}","page":"Models for fluids","title":"MaterialModels.matcalc__density__p","text":"@matcalc(:pressure, model::NewtonianFluid{EOS}; p::Real)\n\nCompute the mass density based on the equation of state EOS.\n\n\n\n\n\n","category":"method"},{"location":"FluidModel/#MaterialModels.matcalc__stress__d__ρ-Tuple{NewtonianFluid, SymmetricSecondOrderTensor{3}, Real}","page":"Models for fluids","title":"MaterialModels.matcalc__stress__d__ρ","text":"@matcalc(:stress, model::NewtonianFluid; d::SymmetricSecondOrderTensor{3}, ρ::Real)\n\nCompute the stress.\n\n\n\n\n\n","category":"method"},{"location":"LinearElastic/#Linear-elastic-model","page":"Linear elastic model","title":"Linear elastic model","text":"","category":"section"},{"location":"LinearElastic/","page":"Linear elastic model","title":"Linear elastic model","text":"LinearElastic","category":"page"},{"location":"LinearElastic/#MaterialModels.LinearElastic","page":"Linear elastic model","title":"MaterialModels.LinearElastic","text":"LinearElastic(; parameters...)\n\nParameters\n\nChoose only 2 parameters.\n\nE: Young's modulus\nK: bulk modulus\nG: shear modulus\nλ: Lamé's first parameter\nν: Poisson's ratio\n\n\n\n\n\n","category":"type"},{"location":"LinearElastic/#Methods","page":"Linear elastic model","title":"Methods","text":"","category":"section"},{"location":"LinearElastic/","page":"Linear elastic model","title":"Linear elastic model","text":"Modules = [MaterialModels]\nOrder = [:function]\nPages = [\"LinearElastic.jl\"]","category":"page"},{"location":"LinearElastic/#MaterialModels.matcalc__stiffness__-Tuple{LinearElastic}","page":"Linear elastic model","title":"MaterialModels.matcalc__stiffness__","text":"@matcalc(:stiffness, model::LinearElastic)\n\nReturn fourth-order stiffness tensor.\n\n\n\n\n\n","category":"method"},{"location":"LinearElastic/#MaterialModels.matcalc__strain__σ-Tuple{LinearElastic, SymmetricSecondOrderTensor{3}}","page":"Linear elastic model","title":"MaterialModels.matcalc__strain__σ","text":"@matcalc(:strain, model::LinearElastic; σ::SymmetricSecondOrderTensor{3})\n\nCompute strain.\n\n\n\n\n\n","category":"method"},{"location":"LinearElastic/#MaterialModels.matcalc__stress__dϵ__σ-Tuple{LinearElastic, SymmetricSecondOrderTensor{3}, SymmetricSecondOrderTensor{3}}","page":"Linear elastic model","title":"MaterialModels.matcalc__stress__dϵ__σ","text":"@matcalc(:stress, model::LinearElastic; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})\n\nCompute stress.\n\n\n\n\n\n","category":"method"},{"location":"LinearElastic/#MaterialModels.matcalc__stress__ϵ-Tuple{LinearElastic, SymmetricSecondOrderTensor{3}}","page":"Linear elastic model","title":"MaterialModels.matcalc__stress__ϵ","text":"@matcalc(:stress, model::LinearElastic; ϵ::SymmetricSecondOrderTensor{3})\n\nCompute stress.\n\n\n\n\n\n","category":"method"},{"location":"SoilHypoelastic/#Hypoelastic-model-for-soils","page":"Hypoelastic model for soils","title":"Hypoelastic model for soils","text":"","category":"section"},{"location":"SoilHypoelastic/","page":"Hypoelastic model for soils","title":"Hypoelastic model for soils","text":"SoilHypoelastic","category":"page"},{"location":"SoilHypoelastic/#MaterialModels.SoilHypoelastic","page":"Hypoelastic model for soils","title":"MaterialModels.SoilHypoelastic","text":"SoilHypoelastic(; κ::Real, ν::Real, e_0::Real)\n\nEquation\n\ndotbmsigma = 3K mathrmvol(dotbmepsilon) + 2G mathrmdev(dotbmepsilon)\n\nwhere\n\nK = frac1+e_0kappa p quad G = frac3(1-nu)2(1+nu)\n\nParameters\n\nκ: elastic compressibility index\nν: Poisson's ratio\ne_0: initial void ratio\n\n\n\n\n\n","category":"type"},{"location":"SoilHypoelastic/#Methods","page":"Hypoelastic model for soils","title":"Methods","text":"","category":"section"},{"location":"SoilHypoelastic/","page":"Hypoelastic model for soils","title":"Hypoelastic model for soils","text":"Modules = [MaterialModels]\nOrder = [:function]\nPages = [\"SoilHypoelastic.jl\"]","category":"page"},{"location":"SoilHypoelastic/#MaterialModels.matcalc__bulkmodulus__σ-Tuple{SoilHypoelastic, SymmetricSecondOrderTensor{3}}","page":"Hypoelastic model for soils","title":"MaterialModels.matcalc__bulkmodulus__σ","text":"@matcalc(:bulkmodulus, model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3})\n\nCompute the bulk modulus.\n\ndotp = K dotepsilon_mathrmv\n\n\n\n\n\n","category":"method"},{"location":"SoilHypoelastic/#MaterialModels.matcalc__shearmodulus__σ-Tuple{SoilHypoelastic, SymmetricSecondOrderTensor{3}}","page":"Hypoelastic model for soils","title":"MaterialModels.matcalc__shearmodulus__σ","text":"@matcalc(:shearmodulus, model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3})\n\nCompute the shear modulus.\n\ndotq = 3G dotepsilon_mathrms\n\n\n\n\n\n","category":"method"},{"location":"SoilHypoelastic/#MaterialModels.matcalc__stiffness__σ-Tuple{SoilHypoelastic, SymmetricSecondOrderTensor{3}}","page":"Hypoelastic model for soils","title":"MaterialModels.matcalc__stiffness__σ","text":"@matcalc(:stiffness, model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3})\n\nCompute the 4th-order stiffness tensor.\n\n\n\n\n\n","category":"method"},{"location":"SoilHypoelastic/#MaterialModels.matcalc__stress__dϵ__σ-Tuple{SoilHypoelastic, SymmetricSecondOrderTensor{3}, SymmetricSecondOrderTensor{3}}","page":"Hypoelastic model for soils","title":"MaterialModels.matcalc__stress__dϵ__σ","text":"@matcalc(:stress, model::SoilHypoelastic; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})\n\nCompute the stress.\n\n\n\n\n\n","category":"method"},{"location":"DruckerPrager/#Drucker-Prager-model","page":"Drucker-Prager model","title":"Drucker-Prager model","text":"","category":"section"},{"location":"DruckerPrager/","page":"Drucker-Prager model","title":"Drucker-Prager model","text":"DruckerPrager","category":"page"},{"location":"DruckerPrager/#MaterialModels.DruckerPrager","page":"Drucker-Prager model","title":"MaterialModels.DruckerPrager","text":"DruckerPrager(::ElasticModel; A, B, b)\n\nYield function\n\nf =  bms  - (A - Bp)\n\nPlastic flow\n\ng =  bms  + b p\n\n\n\n\n\nDruckerPrager(::ElasticModel, mohr_coulomb_type; c, ϕ, ψ = ϕ, tensioncutoff = :auto)\n\nParameters\n\nmohr_coulomb_type: choose from :compression, :tension, :average, :inscribed and :planestrain\nc: cohesion\nϕ: internal friction angle (radian)\nψ: dilatancy angle (radian)\ntensioncutoff: set limit of mean stress, false or :auto (use the mean stress at  bms  = 0)\n\n\n\n\n\n","category":"type"},{"location":"DruckerPrager/#Methods","page":"Drucker-Prager model","title":"Methods","text":"","category":"section"},{"location":"DruckerPrager/","page":"Drucker-Prager model","title":"Drucker-Prager model","text":"Modules = [MaterialModels]\nOrder = [:function]\nPages = [\"DruckerPrager.jl\"]","category":"page"},{"location":"DruckerPrager/#MaterialModels.matcalc__plasticflow__σ-Tuple{DruckerPrager, SymmetricSecondOrderTensor{3}}","page":"Drucker-Prager model","title":"MaterialModels.matcalc__plasticflow__σ","text":"@matcalc(:plasticflow, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})\n\nCompute the plastic flow (the gradient of plastic potential function in stress space, i.e., partialgpartialsigma).\n\n\n\n\n\n","category":"method"},{"location":"DruckerPrager/#MaterialModels.matcalc__stress__dϵ__σ-Tuple{DruckerPrager{LinearElastic}, SymmetricSecondOrderTensor{3}, SymmetricSecondOrderTensor{3}}","page":"Drucker-Prager model","title":"MaterialModels.matcalc__stress__dϵ__σ","text":"@matcalc(:stress, model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})\n\nCompute the stress.\n\n\n\n\n\n","category":"method"},{"location":"DruckerPrager/#MaterialModels.matcalc__stressall__dϵ__σ-Tuple{DruckerPrager{LinearElastic}, SymmetricSecondOrderTensor{3}, SymmetricSecondOrderTensor{3}}","page":"Drucker-Prager model","title":"MaterialModels.matcalc__stressall__dϵ__σ","text":"@matcalc(:stressall, model::DruckerPrager{LinearElastic}; σ::SymmetricSecondOrderTensor{3}, dϵ::SymmetricSecondOrderTensor{3})\n\nCompute the stress and related variables as NamedTuple.\n\n\n\n\n\n","category":"method"},{"location":"DruckerPrager/#MaterialModels.matcalc__yieldfunction__σ-Tuple{DruckerPrager, SymmetricSecondOrderTensor{3}}","page":"Drucker-Prager model","title":"MaterialModels.matcalc__yieldfunction__σ","text":"@matcalc(:yieldfunction, model::DruckerPrager; σ::SymmetricSecondOrderTensor{3})\n\nCompute the yield function.\n\n\n\n\n\n","category":"method"},{"location":"#MaterialModels","page":"Home","title":"MaterialModels","text":"","category":"section"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/KeitaNakamura/MaterialModels.jl.git","category":"page"},{"location":"","page":"Home","title":"Home","text":"or","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> registry add https://github.com/KeitaNakamura/KeitaNakamuraRegistry.git\n\npkg> add MaterialModels","category":"page"},{"location":"utils/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"utils/","page":"Utilities","title":"Utilities","text":"Modules = [MaterialModels]\nOrder = [:function]\nPages = [\"misc.jl\"]","category":"page"},{"location":"utils/#MaterialModels.matcalc__jaumann2caucy__W__dσ_jaumann__σ-Tuple{Tensor{Tuple{3, 3}, T, 2} where T, SymmetricSecondOrderTensor{3}, SymmetricSecondOrderTensor{3}}","page":"Utilities","title":"MaterialModels.matcalc__jaumann2caucy__W__dσ_jaumann__σ","text":"@matcalc(:jaumann2caucy; dσ_jaumann::SymmetricSecondOrderTensor{3}, σ::SymmetricSecondOrderTensor{3}, W::SecondOrderTensor{3})\n\nConvert the Jaumann stress rate to the Caucy stress rate. For example, this function can be used as follows:\n\ndσᴶ = @matcalc(:stress, model; σ, dϵ = symmetric(∇v*dt)) - σ\ndσ = @matcalc(:jaumann2caucy; dσ_jaumann = dσᴶ, σ, W = skew(∇v*dt))\n\nEquation\n\ndotbmsigma^mathrmJ = dotbmsigma - bmW cdot bmsigma + bmsigma cdot bmW\n\n\n\n\n\n","category":"method"},{"location":"utils/#MaterialModels.matcalc__soundspeed__E__ν__ρ-Tuple{Real, Real, Real}","page":"Utilities","title":"MaterialModels.matcalc__soundspeed__E__ν__ρ","text":"@matcalc(:soundspeed; K::Real, G::Real, ρ::Real)\n\nCompute the speed of sound in solids.\n\nEquation\n\nc = sqrtfracE(1-nu)rho (1+nu)(1-2nu)\n\nParameters\n\nE: Young's modulus\nν: Poisson's ratio\nρ: density\n\n\n\n\n\n","category":"method"},{"location":"utils/#MaterialModels.matcalc__soundspeed__G__K__ρ-Tuple{Real, Real, Real}","page":"Utilities","title":"MaterialModels.matcalc__soundspeed__G__K__ρ","text":"@matcalc(:soundspeed; K::Real, G::Real, ρ::Real)\n\nCompute the speed of sound in solids.\n\nEquation\n\nc = sqrtfracK + frac43Grho\n\nParameters\n\nK: bulk modulus\nG: shear modulus\nρ: density\n\n\n\n\n\n","category":"method"}]
}
