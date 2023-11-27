module FourierFloquetHill

#importing necessary libraries and packages
using QuadGK, LinearAlgebra, Dates, Elliptic, PlotlyJS
import Elliptic.Jacobi #Note that Julia's elliptic functions use the m notation (m=k^2)

export FFHM

#this function returns the maximum order of a scalar operator
function max_orders(input_op)
    if occursin("D",input_op) #account for cases with "* D"
        input_op = replace(input_op, "* D" => "*D")

        #determine how many terms the user put in their operator

        count = 0;
        D_indices = [];
        for i=1:length(input_op)
            if input_op[i] == 'D'
                count += 1
                push!(D_indices,i)
            end
        end

        if length(input_op[D_indices[end]:end]) > 2
            count += 1
        end

        #collect array of inputted operator orders, collect array of inputted coefficients

        op_orders = []; #initialize list of input operator orders

        D_indices = []; #initialize list of indices where D is
        op_indices = []; #initialize list of indices where + is (but only plusses between coeffs)

        for i=1:count #iterate through each operator coefficient
            if i == 1 #first coefficient
                nextDInd = findfirst('D',input_op); #find the first "D"

                if typeof(nextDInd) == Nothing #if there is no D (then must be 0th order)
                    push!(op_orders,0)
                else
                    push!(D_indices,nextDInd)
                    nextPlusInd = findnext('+', input_op, nextDInd); #find next plus after last D
                    nextMinusInd = findnext('-', input_op, nextDInd); #find next minus after last D
                    if (typeof(nextPlusInd)==Nothing) || (typeof(nextMinusInd)==Nothing) 
                        if typeof(nextPlusInd)==Nothing #if no next plus, then minus must be start of next coefficient
                            opInd = nextMinusInd;
                            push!(op_indices,opInd)
                        else #otherwise next plus must be start of next coefficient
                            opInd = nextPlusInd;
                            push!(op_indices,opInd)
                        end
                    else #there is both a next plus and a next minus
                        if nextPlusInd < nextMinusInd #if next plus comes sooner, it must be the start of next coeff
                            opInd = nextPlusInd;
                            push!(op_indices,opInd)
                        else #if next minus comes sooner, it must be the start of the next coefficient
                            opInd = nextMinusInd
                            push!(op_indices,opInd)
                        end
                    end
                    op_order = parse(Int,input_op[nextDInd+1:opInd-1]); #extract out order of coefficient (# after D)
                    push!(op_orders,op_order)
                end
            else #coefficients following the first one
                nextDInd = findnext('D',input_op,D_indices[end]+1); #find next D

                if typeof(nextDInd) == Nothing #if no D, then coefficient is order 0
                    push!(op_orders,0)
                else #there are D's
                    push!(D_indices,nextDInd)
                    nextPlusInd = findnext('+',input_op,D_indices[end]); #find next plus after last D
                    nextMinusInd = findnext('-',input_op,D_indices[end]); #find next minus after last D
                    if (typeof(nextPlusInd)==Nothing) || (typeof(nextMinusInd)==Nothing)
                        if typeof(nextPlusInd)==Nothing #if no next plus, then minus must be start of next coefficient
                            opInd = nextMinusInd;
                            push!(op_indices,opInd)
                        else #otherwise next plus must be start of next coefficient
                            opInd = nextPlusInd;
                            push!(op_indices,opInd)
                        end
                    else #there is both a next plus and a next minus
                        if nextPlusInd < nextMinusInd #if next plus comes sooner, it must be the start of next coeff
                            opInd = nextPlusInd;
                            push!(op_indices,opInd)
                        else #if next minus comes sooner, it must be the start of the next coefficient
                            opInd = nextMinusInd
                            push!(op_indices,opInd)
                        end
                    end
                    op_order = parse(Int,input_op[nextDInd+1:opInd-1]); #extract out order of coefficient (# after D)
                    push!(op_orders,op_order)
                end
            end
            i += 1 #move on to next coefficient
        end

    else #no D in operator at all so must be order 0
        op_orders = [0]
    end

    max_order = maximum(op_orders);
    
    return max_order #return maximum order
end

#determine the max order out of all inputs in a vector input
function max_vec_order(vec)
    op_orders= []; #initialize list of orders of operators input in the vector
    for i=1:size(vec)[1] #iterate through operator vector inputs
        for j=1:size(vec)[2]
            op = vec[i,j];
            max_order = max_orders(op); #call max_orders on every input
            push!(op_orders,max_order)
        end
    end
    return maximum(op_orders) #return maximum order
end

#transforms scalar input into an operator the code can read
function scalar_transform(input_op,op_print,max_order=0)
    
    if occursin("D",input_op)
            #account for cases with "* D"
        input_op = replace(input_op, "* D" => "*D")

        #determine how many terms the user put in their operator

        count = 0;
        D_indices = [];
        for i=1:length(input_op)
            if input_op[i] == 'D'
                count += 1
                push!(D_indices,i)
            end
        end

        if length(input_op[D_indices[end]:end]) > 2
            count += 1
        end

        #collect array of inputted operator orders, collect array of inputted coefficients

        op_orders = []; #initialize list of input operator orders
        op_coeffs = []; #initialize list of input operator coefficients

        D_indices = []; #initialize list of indices where D is
        op_indices = []; #initialize list of indices where + is (but only plusses between coeffs)

        for i=1:count #iterate through each coefficient
            if i == 1 #deal with first coefficient
                nextDInd = findfirst('D',input_op); #find the first D in the operator (first coeff ends there)

                if typeof(nextDInd) == Nothing #if no D then only one term of order 0
                    push!(op_orders,0)
                    op_coeff = input_op; #coeff must be entire output
                    push!(op_coeffs,op_coeff)
                else #if D present
                    push!(D_indices,nextDInd)
                    #determine index of operation that separates coefficients as done in max_order function
                    nextPlusInd = findnext('+', input_op, nextDInd);
                    nextMinusInd = findnext('-', input_op, nextDInd);
                    if (typeof(nextPlusInd)==Nothing) || (typeof(nextMinusInd)==Nothing)
                        if typeof(nextPlusInd)==Nothing
                            opInd = nextMinusInd;
                            push!(op_indices,opInd)
                        else
                            opInd = nextPlusInd;
                            push!(op_indices,opInd)
                        end
                    else
                        if nextPlusInd < nextMinusInd
                            opInd = nextPlusInd;
                            push!(op_indices,opInd)
                        else
                            opInd = nextMinusInd
                            push!(op_indices,opInd)
                        end
                    end
                    op_order = parse(Int,input_op[nextDInd+1:opInd-1]); #determine operator order (# after D)
                    push!(op_orders,op_order)
                    op_coeff = input_op[1:nextDInd-1] #coeff is beginning of string up until first D
                    if length(op_coeff) < 3 #deal with case where a "ghost" plus or minus 1 is be coeff
                        if occursin("-",op_coeff) #if -, then must be -1
                            op_coeff = "-1"
                        else #otherwise +1
                            op_coeff = "1"
                        end
                    else #otherwise coeff is beginning up until first D
                        op_coeff = input_op[1:nextDInd-2]
                    end
                    push!(op_coeffs,op_coeff)
                end
            else #deal with the rest of the coeffs
                nextDInd = findnext('D',input_op,D_indices[end]+1); #find the next D after the previous one

                if typeof(nextDInd) == Nothing #if no D then order must be zero
                    push!(op_orders,0)
                    op_coeff = input_op[op_indices[end]:end] #coeff is from previous operation to end
                    if occursin("+",op_coeff[1:2]) #if plus is the previous op, then get rid of that
                        op_coeff = input_op[op_indices[end]+1:end]
                        push!(op_coeffs,op_coeff)
                    else #if minus is previous op, then keep that in
                        push!(op_coeffs,op_coeff)
                    end
                else #if there is a D
                    push!(D_indices,nextDInd)
                    op_coeff = input_op[op_indices[end]:nextDInd-1] #coeff is from previous operation to end
                    if length(op_coeff) < 3 #deal with case where a "ghost" plus or minus 1 is be coeff as above
                        if occursin("+",op_coeff) 
                            op_coeff = "1"
                            push!(op_coeffs,op_coeff)
                        elseif occursin("-",op_coeff)
                            op_coeff = "-1"
                            push!(op_coeffs,op_coeff)
                        end
                    else #otherwise coeff is from previous operatrion up until next D
                        op_coeff = input_op[op_indices[end]:nextDInd-2] 
                        if occursin("+",op_coeff[1:2]) #if plus is the previous op, then get rid of that
                            op_coeff = input_op[op_indices[end]+1:nextDInd-2]
                            push!(op_coeffs,op_coeff)
                        else #if minus is previous op, then keep that in
                            push!(op_coeffs,op_coeff)
                        end
                    end
                    nextPlusInd = findnext('+',input_op,D_indices[end]); #find next plus
                    nextMinusInd = findnext('-',input_op,D_indices[end]); #find next minus
                    if (typeof(nextPlusInd)==Nothing) || (typeof(nextMinusInd)==Nothing) 
                        if typeof(nextPlusInd)==Nothing #if no next plus, next minus must separate coefficients
                            opInd = nextMinusInd;
                            push!(op_indices,opInd)
                        else #if no next minus, next plus must separate coefficients
                            opInd = nextPlusInd;
                            push!(op_indices,opInd)
                        end
                    else #if both a next plus and a next minus
                        if nextPlusInd < nextMinusInd #if next plus comes first, it must separate coefficients
                            opInd = nextPlusInd;
                            push!(op_indices,opInd)
                        else #if next minus comes first, it must separate coefficients
                            opInd = nextMinusInd
                            push!(op_indices,opInd)
                        end
                    end
                    op_order = parse(Int,input_op[nextDInd+1:opInd-1]); #extract out order of coefficient (# after D)
                    push!(op_orders,op_order)
                end
            end
            i += 1 #go to next coefficient
        end
        
        #if the user wants to see how the code interpretted their input operator
        if op_print == true
            operator_string = op_coeffs[1]*"*D"*repr(op_orders[1]);
            for i=2:length(op_coeffs)
                operator_string = operator_string*"+"*op_coeffs[i]*"*D"*repr(op_orders[i])
            end
            display(operator_string)
        end
        
    else #if no D appears in operator at all (must be one coeff of order 0)
        op_coeffs = [input_op] 
        op_orders = [0]
        
        #if the user wants to see how the code interpretted their input operator
        if op_print == true
            operator_string = op_coeffs[1]*"*D"*repr(op_orders[1]);
            for i=2:length(op_coeffs)
                operator_string = operator_string*"+"*op_coeffs[i]*"*D"*repr(op_orders[i])
            end
            display(operator_string)
        end
    end

    if max_order == 0 #if scalar problem, just take max of orders
        max_order = maximum(op_orders);
    end
    
    #convert array of string coeffs to array of function coeffs

    op_func_coeffs = [];
    for i=1:length(op_coeffs)
        try
            func_coeff = eval(Meta.parse("x -> " * op_coeffs[i]));
            push!(op_func_coeffs,func_coeff)
        catch err
            error("Issue with operator input")
        end
    end
    
    #zero pad the rest of the coefficients if they aren't there

    op_funcs = [];

    exist_count = 1;
    for i=max_order:-1:0
        if i in op_orders
            push!(op_funcs,op_func_coeffs[exist_count])
            exist_count += 1;
        else
            push!(op_funcs,x->0)
        end
    end
    
    #reverse the order to go from lowest to highest like my code takes in

    final_op = reverse(op_funcs);
    
    return final_op
end

#this function takes in user input and transforms it into an operator the code can read
function operator_transform(input_op,op_print)
    if typeof(input_op) == String #if operator is scalar
        new_op = scalar_transform(input_op,op_print) #transform the operator to one the Hill's method code can read
    
    else #if the operator is a vector
        ops = []
        max_order = max_vec_order(input_op) #get the maximum order out of the elements in the vector
        for i=1:size(input_op)[1] #transform each operator element in the vector
            for j=1:size(input_op)[2]
                op = scalar_transform(input_op[i,j],op_print,max_order)
                push!(ops,op)
            end
        end
        
        #concatentate the transformed vector operator elements to be in the form that the Hill's method code can read
        op_concat = []
        for k=1:length(ops)
            op = ops[k][1]
            for j=2:length(ops[k])
                op = hcat(op,ops[k][j])
            end
            push!(op_concat,op)
        end

        rows = [];
        op_dim = size(input_op)[1]
        for i=1:op_dim
            row = op_concat[((i-1)*op_dim)+1];
            for j=(((i-1)*op_dim)+2):(i*op_dim)
                row = hcat(row,op_concat[j]);
            end
            push!(rows,row);
        end

        new_op = rows[1];
        for i=2:length(rows)
            new_op = vcat(new_op,rows[i])
        end

    end
    
    return new_op
end

function quad_trap(f,a,b,N) 
    ########################################################
    #function to compute integrals using the trapezoid rule
    #inputs:
    #       f - the function you want to integrate
    #       a - lower bound of integration
    #       b - upper bound of integration
    #       N - how many steps in the integration
    #outputs:
    #       int - value of the definite integral
    ########################################################
    h = (b-a)/N
    int = h * ( f(a) + f(b) ) / 2
    for k=1:N-1
        xk = (b-a) * k/N + a
        int = int + h*f(xk)
    end
    return int
end

function Fcoeffs(f_vec,M,N,P,L)
    ##########################################################################################################
    #function to compute the Fourier coefficients
    #inputs:
    #       f_vec - vector containing coefficients of linear operator for a given block matrix of the problem
    #outputs:
    #       fk_coeffs - matrix of Fourier coefficients for a given block matrix of the problem
    ##########################################################################################################
    fk_coeffs = complex(zeros(Int(M+1),Int(2*(floor((2*N)/P)+1)+1))); #initialize matrix of Fourier coefficients
    for k=0:M #iterate through coefficients of linear operator
        for j=-Int((floor((2*N)/P)+1)):Int((floor((2*N)/P)+1)) #iterate through every Fourier coefficient
                integrand = x-> f_vec[Int(k+1)](x)*exp((-1im*2*pi*j*x)/(L)); #compute integrand
                integral=quad_trap(integrand,-(L)/2,(L)/2,10000) #compute integral
                f_hat = (1/(L))*integral; #compute f_hat
                fk_coeffs[Int(k+1),Int(j+(floor((2*N)/P)+1)+1)] = f_hat; #put coefficient in corresponding matrix entry
        end
    end
    return fk_coeffs
end

#performs Hill's Method
function HillMethod(operator,L,N,D,P,op_dim)
    
    if op_dim == 1
        operator = reverse(operator);
        M = length(operator) - 1; #order of linear operator
        matrix_dim = (2*N)+1; #final matrix will be (2N+1)x(2N+1)

        f_vec = [];
        for i=1:length(operator)
            push!(f_vec,operator[i])
        end
        f_vec = reverse(f_vec)

        #computes the Fourier coefficients
        fk_coeffs = complex(zeros(Int(M+1),Int(2*(floor((2*N)/P)+1)+1))); #initialize matrix of Fourier coefficients
        #iterate through each coefficient of the linear operator
        for k=0:M
            #iterate through all values that need Fourier coefficients computed
            for j=-Int((floor((2*N)/P)+1)):Int((floor((2*N)/P)+1))
                integrand = x -> f_vec[k+1](x)*exp((-1im*2*π*j*x)/(L)); #calculate integrand
                integral=quad_trap(integrand,-(L)/2,(L)/2,10000) #compute integral
                f_hat = (1/(L))*integral; #compute f_hat
                fk_coeffs[k+1,Int(j+(floor((2*N)/P)+1)+1)] = f_hat; #put coefficient in appropriate slot in matrix
            end
        end

        #define range of Floquet parameters
        μ_min = (-π/(P*L));
        μ_max = π/(2*L);
        μ_step = (μ_max-μ_min)/D;

        #initialize arrays for the eigenvalues and eigenvectors
        evals_arr = []
        evecs_arr = []

        #loop through Floquet parameters
        for val = 1:D
            μ = μ_min + (val-1)*μ_step; #compute corresponding Floquet value
            L_matrix = complex(zeros(matrix_dim,matrix_dim)); #initialize (complex) matrix
            #iterate through matrix values
            for m=-N:N
                for n=-N:N
                    #iterate through Fourier coefficients
                    for k=0:M
                        if (mod(n-m,P) == 0) #if n-m even
                            j = Int((n-m)/P); #compute j
                            f_coeff = fk_coeffs[k+1,Int((floor((2*N)/P)+1)+1+j)];  #extract out corresponding Fourier coeff.
                            L_matrix[n+(N+1),m+(N+1)] += f_coeff*((1im*(μ+((2*π*m)/(P*L))))^k) #add element to L matrix
                        end
                    end
                end 
            end

            #compute evals and evecs using eigen
            data = eigen(L_matrix)
            evals = data.values
            push!(evals_arr,evals)
            evecs = data.vectors
            push!(evecs_arr,evecs)
        end

        evals_arr = vcat(evals_arr...)

        #plot spectrum
        display(plot(scatter(x=real(evals_arr),y=imag(evals_arr),mode="markers"),Layout(title="Spectrum")))
        
        return(evals_arr,evecs_arr)

    else
        M = (size(operator)[2]/op_dim) - 1; #order of linear operator
        matrix_dim = (2*N)+1; #final matrix will be (2N+1)x(2N+1)

        op_vec = [];
        for i=1:op_dim
            for j=1:Int((op_dim*(M+1)))
                push!(op_vec,operator[i,j])
            end
        end

        f_mat = Array{Function}(undef, op_dim^2, Int(M+1))
        for i=1:op_dim^2
            f_mat[i,:] = op_vec[((i-1)*(Int(M+1)))+1:(i*(Int(M+1)))]
        end

        Fk_coeffs = [];
        for i=1:op_dim^2
            push!(Fk_coeffs,Fcoeffs(f_mat[i,:],M,N,P,L));
        end

        #define range of Floquet parameters
        μ_min = (-π/(P*L));
        μ_max = π/(2*L);
        μ_step = (μ_max-μ_min)/D;

        #initialize arrays for the eigenvalues and eigenvectors
        evals_arr = []
        evecs_arr = []

        #construct truncated bi-infinite matrix and compute evals/evecs
        for val = 1:D #loop through Floquet parameters
            μ = μ_min + (val-1)*μ_step; #compute corresponding Floquet parameters

            L_blocks = [];
            for i=1:op_dim^2
                L_block = complex(zeros(matrix_dim,matrix_dim)) #initialize (complex) matrix;
                #iterate through matrix values
                for m=-N:N
                    for n=-N:N
                        #iterate through coefficients of linear operators
                        for k=0:M
                            if (mod(n-m,P) == 0) #if n-m even
                                j = Int((n-m)/P); #compute j

                                #extract out corresponding Fourier coefficient
                                f_coeff = Fk_coeffs[i][Int(k+1),Int((floor((2*N)/P)+1)+1+j)];

                                #compute values for each block matrix
                                L_block[n+(N+1),m+(N+1)] += f_coeff*(1im*(μ+((2*π*m)/(P*L))))^k;
                            end
                        end
                    end
                end
                push!(L_blocks,L_block)
            end

            rows = [];
            for i=1:op_dim
                row = L_blocks[((i-1)*op_dim)+1];
                for j=(((i-1)*op_dim)+2):(i*op_dim)
                    row = hcat(row,L_blocks[j]);
                end
                push!(rows,row);
            end

            L_matrix = rows[1];
            for i=2:length(rows)
                L_matrix = vcat(L_matrix,rows[i])
            end

            #compute evals and evecs using eigen
            data = eigen(L_matrix)
            evals = data.values
            push!(evals_arr,evals)
            evecs = data.vectors
            push!(evecs_arr,evecs)
        end

        evals_arr = vcat(evals_arr...)

        #plot the spectrum
        display(plot(scatter(x=real(evals_arr),y=imag(evals_arr),mode="markers"),Layout(title="Spectrum")))

        return(evals_arr,evecs_arr)
        end
    
end

#This function basically does exactly that of the HillMethod code but does this simultaneously on the two operators
# in the generalized EVP. It then computes the generalized evals and evecs using eigen(op1,op2)
function HillMethodGen(operator1,operator2,L,N,D,P,op_dim)
    
    if op_dim == 1
        operator1 = reverse(operator1);
        operator2 = reverse(operator2);
        M1 = length(operator1) - 1; #order of linear operators
        M2 = length(operator2) - 1; #order of second operator
        matrix_dim = (2*N)+1; #final matrices will be (2N+1)x(2N+1)

        f_vec1 = [];
        f_vec2 = [];
        for i=1:length(operator1)
            push!(f_vec1,operator1[i])
        end
        for i=1:length(operator2)
            push!(f_vec2,operator2[i])
        end
        f_vec1 = reverse(f_vec1)
        f_vec2 = reverse(f_vec2)

        #computes the Fourier coefficients
        fk_coeffs1 = complex(zeros(Int(M1+1),Int(2*(floor((2*N)/P)+1)+1))); #initialize matrix of Fourier coefficients
        fk_coeffs2 = complex(zeros(Int(M2+1),Int(2*(floor((2*N)/P)+1)+1)));
        #iterate through each coefficient of the linear operator
        for k=0:M1
            #iterate through all values that need Fourier coefficients computed
            for j=-Int((floor((2*N)/P)+1)):Int((floor((2*N)/P)+1))
                integrand1 = x -> f_vec1[k+1](x)*exp((-1im*2*π*j*x)/(L)); #calculate integrand
                integral1 = quad_trap(integrand1,-(L)/2,(L)/2,10000); #compute integral
                f_hat1 = (1/(L))*integral1; #compute f_hat
                fk_coeffs1[k+1,Int(j+(floor((2*N)/P)+1)+1)] = f_hat1; #put coefficient in appropriate slot in matrix
            end
        end
        for k=0:M2 #repeat for second operator
            #iterate through all values that need Fourier coefficients computed
            for j=-Int((floor((2*N)/P)+1)):Int((floor((2*N)/P)+1))
                integrand2 = x -> f_vec2[k+1](x)*exp((-1im*2*π*j*x)/(L));
                integral2 = quad_trap(integrand2,-(L)/2,(L)/2,10000);
                f_hat2 = (1/(L))*integral2;
                fk_coeffs2[k+1,Int(j+(floor((2*N)/P)+1)+1)] = f_hat2;
            end
        end

        #define range of Floquet parameters
        μ_min = (-π/(P*L));
        μ_max = π/(2*L);
        μ_step = (μ_max-μ_min)/D;

        #initialize arrays for the eigenvalues and eigenvectors
        evals_arr = []
        evecs_arr = []

        #loop through Floquet parameters
        for val = 1:D
            μ = μ_min + (val-1)*μ_step; #compute corresponding Floquet value
            L_matrix = complex(zeros(matrix_dim,matrix_dim)); #initialize (complex) matrix
            T_matrix = complex(zeros(matrix_dim,matrix_dim));
            #iterate through matrix values
            for m=-N:N
                for n=-N:N
                    #iterate through Fourier coefficients
                    for k=0:M1
                        if (mod(n-m,P) == 0) #if n-m even
                            j = Int((n-m)/P); #compute j
                            f_coeff_L = fk_coeffs1[k+1,Int((floor((2*N)/P)+1)+1+j)];  #extract out corresponding Fourier coeff.
                            L_matrix[n+(N+1),m+(N+1)] += f_coeff_L*((1im*(μ+((2*π*m)/(P*L))))^k) #add element to L matrix
                        end
                    end
                    for k=0:M2 #repeat for second operator
                        if (mod(n-m,P) == 0) #if n-m even
                            j = Int((n-m)/P); #compute j
                            f_coeff_T = fk_coeffs2[k+1,Int((floor((2*N)/P)+1)+1+j)];
                            T_matrix[n+(N+1),m+(N+1)] += f_coeff_T*((1im*(μ+((2*π*m)/(P*L))))^k)
                        end
                    end
                end 
            end

            #compute evals and evecs using eigen
            data = eigen(L_matrix,T_matrix)
            evals = data.values
            push!(evals_arr,evals)
            evecs = data.vectors
            push!(evecs_arr,evecs)
        end

        evals_arr = vcat(evals_arr...)

        #plot spectrum
        display(plot(scatter(x=real(evals_arr),y=imag(evals_arr),mode="markers"),Layout(title="Spectrum")))
        
        return(evals_arr,evecs_arr)

    else
        M1 = (size(operator1)[2]/op_dim) - 1; #order of linear operator
        M2 = (size(operator2)[2]/op_dim) - 1; #order of second operator
        matrix_dim = (2*N)+1; #final matrix will be (2N+1)x(2N+1)

        op1_vec = [];
        op2_vec = [];
        for i=1:op_dim
            for j=1:Int((op_dim*(M1+1)))
                push!(op1_vec,operator1[i,j])
            end
            for j=1:Int((op_dim*(M2+1))) #repeat for second operator
                push!(op2_vec,operator2[i,j])
            end
        end

        f1_mat = Array{Function}(undef, op_dim^2, Int(M1+1))
        f2_mat = Array{Function}(undef, op_dim^2, Int(M2+1))
        for i=1:op_dim^2
            f1_mat[i,:] = op1_vec[((i-1)*(Int(M1+1)))+1:(i*(Int(M1+1)))]
            f2_mat[i,:] = op2_vec[((i-1)*(Int(M2+1)))+1:(i*(Int(M2+1)))]
        end

        Fk1_coeffs = [];
        Fk2_coeffs = [];
        for i=1:op_dim^2
            push!(Fk1_coeffs,Fcoeffs(f1_mat[i,:],M1,N,P,L));
            push!(Fk2_coeffs,Fcoeffs(f2_mat[i,:],M2,N,P,L));
        end

        #define range of Floquet parameters
        μ_min = (-π/(P*L));
        μ_max = π/(2*L);
        μ_step = (μ_max-μ_min)/D;

        #initialize arrays for the eigenvalues and eigenvectors
        evals_arr = []
        evecs_arr = []

        #construct truncated bi-infinite matrix and compute evals/evecs
        for val = 1:D #loop through Floquet parameters
            μ = μ_min + (val-1)*μ_step; #compute corresponding Floquet parameters

            L_blocks = [];
            T_blocks = [];
            for i=1:op_dim^2
                L_block = complex(zeros(matrix_dim,matrix_dim)) #initialize (complex) matrix;
                T_block = complex(zeros(matrix_dim,matrix_dim))
                #iterate through matrix values
                for m=-N:N
                    for n=-N:N
                        #iterate through coefficients of linear operators
                        for k=0:M1
                            if (mod(n-m,P) == 0) #if n-m even
                                j = Int((n-m)/P); #compute j

                                #extract out corresponding Fourier coefficient
                                f1_coeff = Fk1_coeffs[i][Int(k+1),Int((floor((2*N)/P)+1)+1+j)];

                                #compute values for each block matrix
                                L_block[n+(N+1),m+(N+1)] += f1_coeff*(1im*(μ+((2*π*m)/(P*L))))^k;
                            end
                        end
                        for k=0:M2 #repeat for second operator
                            if (mod(n-m,P) == 0) #if n-m even
                                j = Int((n-m)/P); #compute j

                                #extract out corresponding Fourier coefficient
                                f2_coeff = Fk2_coeffs[i][Int(k+1),Int((floor((2*N)/P)+1)+1+j)];

                                #compute values for each block matrix
                                T_block[n+(N+1),m+(N+1)] += f2_coeff*(1im*(μ+((2*π*m)/(P*L))))^k;
                            end
                        end
                    end
                end
                push!(L_blocks,L_block)
                push!(T_blocks,T_block)
            end

            rows_L = [];
            rows_T = [];
            for i=1:op_dim
                row_L = L_blocks[((i-1)*op_dim)+1];
                row_T = T_blocks[((i-1)*op_dim)+1];
                for j=(((i-1)*op_dim)+2):(i*op_dim)
                    row_L = hcat(row_L,L_blocks[j]);
                    row_T = hcat(row_T,T_blocks[j]);
                end
                push!(rows_L,row_L);
                push!(rows_T,row_T);
            end

            L_matrix = rows_L[1];
            T_matrix = rows_T[1];
            for i=2:length(rows_L)
                L_matrix = vcat(L_matrix,rows_L[i])
                T_matrix = vcat(T_matrix,rows_T[i])
            end

            #compute evals and evecs using eigen
            data = eigen(L_matrix,T_matrix)
            evals = data.values
            push!(evals_arr,evals)
            evecs = data.vectors
            push!(evecs_arr,evecs)
        end

        evals_arr = vcat(evals_arr...)

        #plot the spectrum
        display(plot(scatter(x=real(evals_arr),y=imag(evals_arr),mode="markers"),Layout(title="Spectrum")))

        return(evals_arr,evecs_arr)
        end
    
end

#takes in user input operator(s) and corresponding info and returns evals and evecs with interactive plot
function FFHM(scalar_period,scalar_num_FourierModes,scalar_num_FloqVals,P,scalar_op_dim,op_print,input_op1,input_op2=0)
    if input_op2 == 0 #if there is only one operator
        new_op = operator_transform(input_op1,op_print); #transform user input operator to one Hill's method code can read
        #call Hill's method code
        (evals_arr,evecs_arr) = Base.invokelatest(HillMethod,new_op,scalar_period,scalar_num_FourierModes,scalar_num_FloqVals,P,scalar_op_dim)
        return (evals_arr,evecs_arr)
    else #generalized EVP (two operators)
        new_op1 = operator_transform(input_op1,op_print); #transform user input operator to one Hill's method code can read
        new_op2 = operator_transform(input_op2,op_print); #transform user input operator to one Hill's method code can read
        #call Hill's method code
        (evals_arr,evecs_arr) = Base.invokelatest(HillMethodGen,new_op1,new_op2,scalar_period,scalar_num_FourierModes,scalar_num_FloqVals,P,scalar_op_dim)
        return (evals_arr,evecs_arr)
    end
end

end #module


