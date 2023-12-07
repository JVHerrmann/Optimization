#  TRABALHO 2 DE OTIMIZAÇÃO ESTRUTURAL - UDESC - PPG Eng. Mecânica - 2023 / 2
#  João Vítor Herrmann
#
# PROBLEMA: Placa plana submetida a tração com furo inserido. 
#           Algumas variáveis de projeto são apresentadas e 
#
# SITUAÇÃO 1 - APENAS r do furo como variável de projeto
#
# Min c´x
#
# A x <= b
#
#
# c = \nabla f(xk)
#
# A = [ \nabla g_j (xk) ' ]
#
# b = \nabla g_j(xk)' * xk - g_j(xk)
#
#


using JuMP
using LinearAlgebra
using LinearAlgebraX
using Plots
import Tulip

# Método de Otimização ser utilizado: SLP - Abaixo o código do método:
function LP_profissional(c::Vector,A::Matrix,b::Vector,ci::Vector,cs::Vector,verbose=false)

    # Dimensões do problema
    n = length(c)
    m = length(b)
    
    size(A)==(m,n) || throw("LP::A deve ser $m × $n")
    length(ci)==length(cs)==n || throw("LP::ci e cs devem ter dimensão $n ")

    # Inicializa o modelo
    lp = Model(Tulip.Optimizer)

    # Define as variáveis
    @variable(lp, x[1:n])

    # Impoe as restrições laterais
    set_lower_bound.(x, ci)
    set_upper_bound.(x, cs)

    # Define o objetivo
    @objective(lp, Min, c'*x)

    # Define as restrições 
    @constraint(lp, restricoes, A*x .<= b)
   
    # Se verbose, mostra o problema que está sendo solucionado
    if verbose
       println("Objetivo ",c'*x)
       for i in restricoes
           println(i)
       end
    end

    # Soluciona
    JuMP.optimize!(lp)

    # Verifica o flag de solução
    st = termination_status(lp)
    println("Flag de solução: $st")

    # Acessa a solução
    objetivo = objective_value(lp)
    sol = value.(x)

    return sol, objetivo

end


const t = 5 #mm
const W = 5*t #25mm
const L = 4*W #100mm
const ntracao = 2 
const ub = 0.01*L #
#const utot = 0.01*L

const F = 4500 # N
const E = 169E9 #GPa - Ferro Fundido Cinzento
const sesc = 214E6 #MPa - Ferro Fundido Cinzento


#Considerando fojbetivo sendo o volume

function objetivo(x::Vector)
    r = x[1]
    return t*(L*W-π*r^2)  #minimizando o volume V(r). Volume total = volume do paralelepípedo - Volume do furo
end

function dobjetivo(x::Vector)
    r = x[1]
    dr = -2*π*r*t
    return [dr]
end


function g1(x::Vector)
    r = x[1]
    -ub+(F/(E*t))*((L-2*r)/W + 2*r/(W-2*r)) #restrição quanto deslocamento máximo ub
end


function dg1(x::Vector)
    r = x[1]
    dr = -(8*F*r*(r-W))/(E*t*W*(W-2*r)^2)
    [dr]
end

function g2(x::Vector)
    r = x[1]
    (sesc/((3-3.13*(2*r/W)+3.66*(2*r/W)^2-1.53*(2*r/W)^3)*F/(2*t*(W-2*r))))-ntracao #restrição quanto ao coeficiente de segurança para a vida infinita de fadiga. Considerando o Kt(r) com aproximação polinomial. 
end

function dg2(x::Vector)
    r = x[1]
    dr = (sesc*t*W^3*(-0.653595*r^3 + 0.881071*r^2*W - 0.390875*r*W^2 + 0.00347089*W^3))/(F*(r^3 - 1.19608*r^2*W + 0.511438*r*W^2 - 0.245098*W^3)^2)
    [dr]
end

function Lineariza(xk::Vector, lm=1E-3)


   # O vetor c é igual a derivada do objetivo no ponto xk
   c = dobjetivo(xk)

   # A matriz A tem as derivadas das restrições em cada em cada linha
   A = [ dg1(xk)' ;
         dg2(xk)' ]

   b = [ dot(dg1(xk),xk) - g1(xk) ;
         dot(dg2(xk),xk) - g2(xk)]    

   # Restrições laterais originais do problema
   
   #xmin = 1E-2*ones(2)  
   xmin = 1*ones(1)  
    
   #xmax = 1E-1*ones(2)  - Valor usado no problema da viga
   xmax = W/2*ones(1) # r do furo não pode ser maior do que metade da largura 

   # Limites móveis, evitando violar as restrições 
   # laterais originais
   ci = max.(xmin,min.(xmax,xk.-lm))
   cs = max.(xmin,min.(xmax,xk.+lm))

   return c,A,b,ci,cs

end




function SLP_Placa_com_furo()
    
    # Define um ponto inicial para começarmos o SLP
    x1 = [2]
 
    # seta um valor de limite móvel
    lm = 1E-2

    # Iterações do SLP
    for i=1:100
 
     # Lineariza o problema no entorno de x1. 
     c,A,b,ci,cs = Lineariza(x1,lm)
 
     # Teste interno para ver se os limites móveis fazem sentido
     @assert all(ci.<=x1.<=cs) 
 
     # Soluciona o problema linearizado
     xopt, fopt = LP_profissional(c,A,b,ci,cs,false)      
 
    # Verifica se a solução do LP satisfaz os limites móveis.
    # Caso contrário, coloca o limite móvel.
     for i in LinearIndices(xopt)
         if xopt[i]<ci[i]
             xopt[i] = ci[i]
         end
         if xopt[i]>cs[i]
             xopt[i] = cs[i]
         end
         
     end
     
    # Calcula a norma da variação do x
    delta = norm(xopt.-x1)

    println("Variação da norma do x foi ", delta)
    @show g1(xopt)
    @show g2(xopt)

    # Copia para a solução atual 
    x1 .= xopt
    @show x1
    
    
    end
 
    return x1
    
 end


SLP_Placa_com_furo()


