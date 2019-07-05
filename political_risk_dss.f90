!Replicate the anaysis of Cuadra & Sapriza 2008
module param
DOUBLE PRECISION :: beta, sigma, r, b_inf, b_sup, y_inf, y_sup, eps_inf, eps_sup, std_eps,&
                    rho, mean_y, std_y, pi_number, lambda, zero, width, prob_excl_end, mean_g,&
                    scale_y
!PARAMETERIZATION IN AGUIAR AND GOPINATH

parameter (beta = 0.953d+0, sigma = 2d+0, r = 0.017d+0, &
           std_eps = 0.025d+0, rho = 0.945d+0, mean_y = 0d+0,  pi_number = 3.1415926535897932d+0, prob_office = 0.90d+0, &
          lambda = 0.969d+0, zero = 0d+0,  width = 1.5d+0, prob_excl_end = 0.282d+0, mean_g = 1.000d+0,&
          scale_y = 10d+0, theta = 0.62d+0, share = (theta/(1d+0 - theta))**(1d+0/sigma))


!y = log(income)
!sigma = coefficient of relative risk aversion
!r = international interest rate
!rho = coefficient of autocorrelation in income
!std_eps = standard deviation of innovation to income shocks
!mean_y = unconditional mean of the growth trend shock
!std_y = standard deviation of growth trend shock

INTEGER :: b_num, y_num, eps_num, quad_num, nout, i_y_global, i_excl_global, i_ss_initial,&
           i_b_global, i_default_global, indicator_2, i_b_next, i_b_zero, b_num_interm
parameter (b_num= 200, y_num =200, b_num_interm=100)

DOUBLE PRECISION :: b_grid(1:b_num), y_grid(1:y_num), default_grid(1:2), indicator_tirar, &
                    matrix_v_global(b_num, y_num), vector_v_global(y_num), y_initial, b_initial, b_global, w1_matrix(y_num),&
                    v1_matrix(y_num), exp_v_vector(b_num), q_vector(b_num), exp_v_excl_vector(b_num), exp_w_vector(b_num), exp_w_excl_vector(b_num)

!y_initial = growth rate state. USED TO SOLVE OPTIMAL SAVINGS RULE AT EVERY ITERATION
!b_initial = debt state. USED TO SOLVE OPTIMAL SAVINGS RULE AT EVERY ITERATION

DOUBLE PRECISION, DIMENSION(b_num, y_num) :: v_matrix, w_matrix
INTEGER, DIMENSION(b_num, y_num) :: default_decision, b_next_matrix, b0_next_matrix_indice
DOUBLE PRECISION, DIMENSION(b_num, y_num) :: v0_matrix, b0_next_matrix, b1_next_matrix, w0_matrix, q_matrix_nodef
DOUBLE PRECISION, DIMENSION(y_num, y_num) :: trans_matrix, cdf_matrix


!2nd component = h (default decision in previous period
!4th component = b (outstanding debt)
!5th component = y (current income realization)

end module


!DEFINE UTILITY FUNCTION
DOUBLE PRECISION function u_fun(c)
USE param
DOUBLE PRECISION, INTENT(IN) :: c
DOUBLE PRECISION :: c_min, cons1, cons2

c_min = 1d-10
cons2 = MAX((c/(1d+0+share)),c_min)
cons1 = MAX(share*cons2,c_min)

u_fun = theta*(cons1**(1-sigma)) / (1-sigma) + (1-theta)*(cons2**(1-sigma) ) / (1-sigma)
!WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)')  c, cons1, cons2, u_fun
end

DOUBLE PRECISION function u_fun_w(c)
USE param
DOUBLE PRECISION, INTENT(IN) :: c
DOUBLE PRECISION :: c_min, cons1, cons2

c_min = 1d-10
cons1 = MAX((c/(1d+0+share)),c_min)
cons2 = MAX(share*cons1,c_min)

u_fun_w = theta*(cons1**(1-sigma)) / (1-sigma) + (1-theta)*(cons2**(1-sigma) ) / (1-sigma)
!WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)')  c, cons1, cons2, u_fun_w
end

!SPECIFY GRID VALUES
subroutine compute_grid
USE param
INTEGER :: n
DOUBLE PRECISION :: dos, DNORDF, delta_y, prob_mass, y_left, y_right, b_inf_interm, b_sup_interm
EXTERNAL :: DNORDF

std_y = std_eps/ SQRT(1 - rho**2)

y_inf = mean_y - 4*std_y
y_sup = mean_y + 4*std_y
delta_y = 4*std_y / (y_num - 1d+0)


b_inf = -3.30d+0  !0.65*EXP(y_inf)    !NEED A SYMMETRIC GRID, SO THAT
b_sup =  1.50d+0 !*EXP(y_inf)

open (10, FILE='graphs\b_grid_dss.txt',STATUS='replace')
open (11, FILE='graphs\y_grid_dss.txt',STATUS='replace')
open (13, FILE='graphs\bounds_dss.txt',STATUS='replace')



WRITE(13, '(F15.11, X, F15.11)') b_inf, b_sup
WRITE(13, '(F15.11, X, F15.11)') y_inf, y_sup

b_inf_interm = -0.04
b_sup_interm = 0.0d+0


n = INT((b_num-b_num_interm)*(b_inf_interm - b_inf)/(b_sup - b_inf))
do i=1,n
   b_grid(i) = b_inf + (b_inf_interm - b_inf)*(i-1d+0)/(n - 1d+0)
end do

n = b_num_interm
do i=1,n
   b_grid(INT((b_num-b_num_interm)*(b_inf_interm - b_inf)/(b_sup - b_inf)) + i) = &
   b_inf_interm +.00001d+0 + (b_sup_interm - b_inf_interm-.00002d+0)*(i-1d+0)/(n-1d+0)
end do

i_b_zero = INT((b_num-b_num_interm)*(b_inf_interm - b_inf)/(b_sup - b_inf)) + b_num_interm + 1

n = b_num - INT((b_num-b_num_interm)*(b_inf_interm - b_inf)/(b_sup - b_inf)) - b_num_interm
do i=1,n
   b_grid(INT((b_num-b_num_interm)*(b_inf_interm - b_inf)/(b_sup - b_inf)) + b_num_interm+ i) = &
   b_sup_interm + (b_sup - b_sup_interm)*(i-1d+0)/(n-1d+0)
end do


do i=1,b_num
 WRITE(10, '(F12.8)') b_grid(i)
end do 

do i=1,y_num
   y_grid(i) = y_inf + (y_sup - y_inf) * (i-1) / (y_num - 1)
   WRITE(11, '(F12.8)') y_grid(i)
end do



!TRANSTION MATRIX AS IN ARELLANO
do i=1,y_num
    y_left  = (y_grid(1) + delta_y - rho*y_grid(i))/ std_eps
    y_right = (y_grid(y_num) - delta_y - rho*y_grid(i))/ std_eps

    trans_matrix(i,1) = DNORDF(y_left)
    trans_matrix(i,y_num) = 1d+0 - DNORDF(y_right)
    !WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)')  y_left, y_right, trans_matrix(i,1)
    do j=2,y_num-1
        y_right = (y_grid(j) + delta_y - rho*y_grid(i))/ std_eps
        y_left  = (y_grid(j) - delta_y - rho*y_grid(i))/ std_eps
        trans_matrix(i,j) = (DNORDF(y_right) - DNORDF(y_left) ) !/prob_may_num
        !WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)')  y_left, y_right, trans_matrix(i,j)
    end do
    !WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)')  y_left, y_right, trans_matrix(i,y_num)
end do

do i=1,y_num
   cdf_matrix(i,1) = trans_matrix(i,1)
   do j=2,y_num
      cdf_matrix(i,j) = cdf_matrix(i,j-1) + trans_matrix(i,j)
   end do
end do

dos = 2d+00

CLOSE(10)
CLOSE(11)
!CLOSE(12)
CLOSE(13)
default_grid(1) = 0d+00
default_grid(2) = 1d+00

end subroutine




DOUBLE PRECISION function q_fun(i_b, i_y)
USE param
INTEGER, INTENT(IN) :: i_b, i_y
DOUBLE PRECISION :: partial_integral, prob_default_same_type
INTEGER :: i, index_tirar

prob_default_other_type  = 0d+0
prob_default_same_type  = 0d+0

do i=1,y_num
    prob_default_same_type  = prob_default_same_type  + default_grid(default_decision(i_b, i)) * trans_matrix(i_y, i)
end do
partial_integral =  1d+0 - prob_default_same_type
q_fun = partial_integral / (1d+0 + r)

end

!Compute the objective function in the Belman equation
DOUBLE PRECISION function objective_function(b)
USE param
DOUBLE PRECISION, INTENT(IN) :: b
INTEGER :: other_type, i,j,h, t, d
DOUBLE PRECISION :: u_fun, q_fun, acum, exp_v_next, exp_v_excl, value_next, y_next, &
                    y_threshold, y_max, y_min, y_next1, acum2,acum1, g, q, DNORDF,&
                    value_next_tirar, acum_int, acum_int_lo, acum_int_hi, output

EXTERNAL u_fun, q_fun, DNORDF



exp_v_next = exp_v_vector(i_b_next)      !0d+0
exp_v_excl = exp_v_excl_vector(i_b_next) !0d+0

d  = i_default_global
q = q_vector(i_b_next)* (1-default_grid(i_default_global)) !q_fun(i_b_next, i_y_global) * mean_g
if (default_grid(i_excl_global) > 0.5 ) then
    output = scale_y * MIN(lambda, exp(y_initial))
else
    output = scale_y * exp(y_initial)
end if

objective_function = u_fun(output + &
                              b_initial * (1-default_grid(i_default_global))- b*q) + beta * mean_g**(1-sigma) * ( &
default_grid(i_excl_global) * (prob_excl_end *exp_v_next + (1-prob_excl_end) * exp_v_excl) + &
(1- default_grid(i_excl_global)) * exp_v_next)

end

!Compute the objective function in the Belman equation
DOUBLE PRECISION function objective_function_w(b)
USE param
DOUBLE PRECISION, INTENT(IN) :: b
INTEGER :: other_type, i,j,h, t, d
DOUBLE PRECISION :: u_fun, q_fun, acum, exp_v_next, exp_v_excl, value_next, y_next, &
                    y_threshold, y_max, y_min, y_next1, acum2,acum1, g, q, DNORDF,&
                    value_next_tirar, acum_int, acum_int_lo, acum_int_hi, output, u_fun_w

EXTERNAL u_fun_w, q_fun, DNORDF


exp_w_next = exp_w_vector(i_b_next)      !0d+0
exp_w_excl = exp_w_excl_vector(i_b_next) !0d+0

d  = i_default_global
q = q_vector(i_b_next)* (1-default_grid(i_default_global)) 

if (default_grid(i_excl_global) > 0.5 ) then
    output = scale_y * MIN(lambda, exp(y_initial))
else
    output = scale_y * exp(y_initial)
end if

objective_function_w = u_fun_w(output + &
                              b_initial * (1-default_grid(i_default_global))- b*q) + beta * ( &
default_grid(i_excl_global) * (prob_excl_end *exp_w_next + (1-prob_excl_end) * exp_w_excl) + &
(1- default_grid(i_excl_global)) * exp_w_next)

end
    
subroutine optimize(index_opt, v_value)
USE param
integer :: MAXFN, t, d, index_opt, index_vector(1)
DOUBLE PRECISION :: b_next, v_value, objective_function, new_value, q_fun, old_value, g, b, u_fun,q, vector(b_num)
EXTERNAL objective_function, q_fun, u_fun

old_value = -10d+3
index_opt = 1
do i=1,b_num
    i_b_next = i
    b_next   = b_grid(i)
    vector(i)  = objective_function(b_next)
    !new_value = objective_function(b_next)
    !if (new_value >= old_value) then
    !    index_opt = i
    !    old_value = new_value
    !end if
end do
!v_value = old_value
index_vector = maxloc(vector)
v_value = vector(index_vector(1))
index_opt = index_vector(1)

10 end subroutine


subroutine iterate
USE param
INTEGER :: d
DOUBLE PRECISION :: y_valor, b_valor, b_valor_def, v_valor, convergence, criteria, deviation, q_fun, b_next, g,&
                    b0_next, b1_next, b, w_value, v_valor_exl, objective_function,  objective_function_w, dev_v, dev_q_paid,&
                    v0_value, v1_value, q, v1_matrix_new(y_num), exp_v_excl, exp_v_next, exp_w_excl, exp_w_next, w1_matrix_new(y_num)

DOUBLE PRECISION, DIMENSION(b_num, y_num) :: v0_matrix_new, w0_matrix_new, q_matrix_nodef_new
DOUBLE PRECISION, DIMENSION(b_num, y_num) :: v_matrix_new, dev_matrix_v, w_matrix_new, dev_matrix_q
INTEGER i_b, i_y, i_def_opt, i_b_optimal0, i_b_optimal1, i_b_optimal
EXTERNAL q_fun, objective_function, objective_function_w




criteria = 1d-5   !CRITERIA FOR CONVERGENCE
convergence = -1
ERRREL = 1d-10    !PRECISION WITH WHICH THE CHEBYCHEV COEFFICIENTS ARE COMPUTED


do WHILE(convergence<0)
    deviation = 0
    dev_v = 0d+0
    dev_q_paid   = 0d+0
    do i_y = 1,y_num
        i_y_global = i_y
        y_initial = y_grid(i_y)
        do i_b = 1,b_num
            exp_v_excl = 0d+0
            exp_v_next = 0d+0
            exp_w_excl = 0d+0
            exp_w_next = 0d+0
            do i=1,y_num
                exp_v_next = exp_v_next + prob_office*v_matrix(i_b, i) * trans_matrix(i_y_global, i) + (1-prob_office)*w_matrix(i_b, i) * trans_matrix(i_y_global, i)
                exp_w_next = exp_w_next + prob_office*w_matrix(i_b, i) * trans_matrix(i_y_global, i) + (1-prob_office)*v_matrix(i_b, i) * trans_matrix(i_y_global, i)
                exp_v_excl = exp_v_excl + prob_office*v1_matrix(i) * trans_matrix(i_y_global, i) + (1-prob_office)*w1_matrix(i) * trans_matrix(i_y_global, i) 
                exp_w_excl = exp_w_excl + prob_office*w1_matrix(i) * trans_matrix(i_y_global, i) + (1-prob_office)*v1_matrix(i) * trans_matrix(i_y_global, i) 
            end do
            exp_v_vector(i_b) = exp_v_next
            exp_v_excl_vector(i_b) = exp_v_excl
            exp_w_vector(i_b) = exp_w_next
            exp_w_excl_vector(i_b) = exp_w_excl
            !WRITE(nout,'(I3, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') i_b, exp_v_next, exp_w_next, exp_v_excl, exp_w_excl
            q_vector(i_b) = q_fun(i_b, i_y_global)
        end do
        i_default_global=2 !Country defaults
        i_excl_global=2
        !WRITE(nout,*)  i_y, i_b_zero
        i_b_next = i_b_zero
        i_b_optimal1 = i_b_zero
        b_next   = b_grid(i_b_next)
        v1_matrix_new(i_y) = objective_function(b_next)
        w1_matrix_new(i_y) = objective_function_w(b_grid(i_b_optimal1))
        
        do i_b = 1, b_num
            i_b_global = i_b
            b_initial = b_grid(i_b)
            i_default_global=1   ! Country does not default and is not excluded for the next period
            i_excl_global=1
            call optimize(i_b_optimal0, v_valor)
            indicator_tirar=0
            b0_next_matrix(i_b, i_y) = b_grid(i_b_optimal0)
            b0_next_matrix_indice(i_b, i_y) = i_b_optimal0
            i_b_next = i_b_optimal0
            v0_matrix_new(i_b, i_y) = v_valor
            w0_matrix_new(i_b, i_y) = objective_function_w(b_grid(i_b_optimal0))
            q_matrix_nodef_new(i_b, i_y) = q_fun(i_b_optimal0, i_y_global)
                   
            !WRITE(nout,'(I3, X, I3, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') i_b, i_b_optimal0, v0_matrix_new(i_b, i_y), w0_matrix_new(i_b, i_y)
            indicator_tirar =0
            if (v1_matrix_new(i_y) <= v0_matrix_new( i_b, i_y)) then
                default_decision( i_b, i_y) = 1
                b_next_matrix(i_b, i_y) = i_b_optimal0  !SAVINGS IF NOT EXCLUDED
                v_matrix_new(i_b, i_y) = v0_matrix_new( i_b, i_y)
                w_matrix_new(i_b, i_y) = w0_matrix_new( i_b, i_y)
            else
                default_decision(i_b, i_y) = 2
                b_next_matrix(i_b, i_y) = i_b_optimal1
                v_matrix_new(i_b, i_y) = v1_matrix_new(i_y)
                w_matrix_new(i_b, i_y) = w1_matrix_new(i_y)
            end if
            dev_q_paid =   MAX(ABS(q_matrix_nodef_new(i_b, i_y) - q_matrix_nodef(i_b, i_y)), dev_q_paid)
            dev_v = max(ABS(v_matrix_new( i_b, i_y) - v_matrix(i_b, i_y)), dev_v)
            deviation = MAX(deviation, MAX(dev_v, dev_q_paid))
            !deviation = MAX(deviation, dev_v)
            dev_matrix_v(i_b, i_y) = ABS(v_matrix_new(i_b, i_y) - v_matrix( i_b, i_y))
            dev_matrix_q(i_b, i_y) = ABS(q_matrix_nodef_new(i_b, i_y) - q_matrix_nodef(i_b, i_y))
            !WRITE(nout, '(I3, X, I4, X, F12.8)') i_y, i_b, dev_matrix(i_b, i_y)
            indicator_tirar=0
        end do
    end do
    WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') deviation, dev_q_paid, dev_v

!3) SAVE RESULTS OF THE CURRENT ITERATION

      open (10, FILE='graphs\v_dss.txt',STATUS='replace')
      open (11, FILE='graphs\default_dss.txt',STATUS='replace')
      open (12, FILE='graphs\q_dss.txt',STATUS='replace')
      open (13, FILE='graphs\b_next_dss.txt',STATUS='replace')
      open (14, FILE='graphs\dev_dss.txt',STATUS='replace')
      open (16, FILE='graphs\q_paid_dss.txt',STATUS='replace')
      open (17, FILE='graphs\b0_next_dss.txt',STATUS='replace')
      open (18, FILE='graphs\w_dss.txt',STATUS='replace')

    do i_b = 1,b_num
        do i_y = 1,y_num
            i_y_global = i_y
            y_initial = y_grid(i_y)
            b_initial = b_grid(i_b)
            i_b_global = i_b

            indicator_tirar=0
            WRITE(10, '(F15.10, X, F15.10, X, F15.10)') v_matrix_new(i_b, i_y), &
            v0_matrix_new(i_b, i_y), v1_matrix_new( i_y)
            WRITE(11, '(F6.2)') default_grid(default_decision(i_b, i_y))
            d = default_decision(i_b, i_y)
            b = b_grid(i_b)
            g = y_grid(i_y)
            WRITE(12, '(F15.11)') q_fun(i_b, i_y)
            WRITE(13, '(F15.11)') b_grid(b_next_matrix(i_b, i_y))
            WRITE(17, '(F15.11)') b0_next_matrix(i_b, i_y)
            WRITE(14, '(F15.11, X,F15.11)') dev_matrix_v(i_b, i_y), dev_matrix_q(i_b, i_y)

            WRITE(16, '(F15.11)') q_fun(b_next_matrix(i_b, i_y), i_y)
            indicator_tirar=0
            WRITE(18, '(F15.10, X, F15.10, X, F15.10)') w_matrix_new(i_b, i_y), &
            w0_matrix_new(i_b, i_y), w1_matrix_new( i_y)
        end do
    end do
    CLOSE(10)
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
    CLOSE(14)
    CLOSE(16)
    CLOSE(17)
    CLOSE(18)
    !UPDATE VALUES OF MATRICES
    v0_matrix = v0_matrix_new
    v1_matrix = v1_matrix_new
    v_matrix = v_matrix_new
    w0_matrix = w0_matrix_new
    w1_matrix = w1_matrix_new
    w_matrix = w_matrix_new
    q_matrix_nodef = q_matrix_nodef_new
   
    !PRINT*,'deviation =', deviation
    if (deviation < criteria) then
        convergence =1
    end if
end do
!FINALLY, STORE VALUE FUNCTIONS, POLICY FUNCTIONS AND PRICES USING A FINER GRID.
!THIS HELPS TO VISUALIZE THE RESULTS.
end subroutine


program main
include 'link_fnl_static.h'
USE param
DOUBLE PRECISION :: y_valor, b_valor, f_valor, q_fun, u_fun, start_time, end_time, indicator_external, def, u_fun_w
INTEGER  i_b, i_y
EXTERNAL q_fun, u_fun, u_fun_w

call cpu_time(start_time)
call compute_grid


indicator_external = 1 !FROM EXTERNAL FILE
if (indicator_external < 0.5) then
    default_decision = 2
    do i_y = 1,y_num
        y_initial = y_grid(i_y)
        !WRITE(nout, '(A10, X, A10, X, A10, X, A10)') 'y', 'b', 'c1', 'c0'
        v1_matrix( i_y) = u_fun( scale_y *MIN( lambda, exp(y_initial)) )
        w1_matrix( i_y) = u_fun_w( scale_y *MIN( lambda, exp(y_initial)) )
        do i_b = 1,b_num
            b_initial = b_grid(i_b)
            v0_matrix(i_b, i_y) = u_fun(EXP(y_grid(i_y)) + b_grid(i_b))
            w0_matrix(i_b, i_y) = u_fun_w(EXP(y_grid(i_y)) + b_grid(i_b))
            v_matrix(i_b, i_y) = v1_matrix(i_y)
            w_matrix(i_b, i_y) = w1_matrix(i_y)
            !WRITE(nout, '(I3, X, I3, X, F12.8)') i_y, i_b, v0_matrix(i_b, i_y)
            !WRITE(nout, '(F10.6, X, F10.6, X, F10.6, X, F10.6)') EXP(y_grid(i_y)), b_grid(i_b),&
        end do
        !pause
    end do
    default_decision = 2   !IN THE LAST PERIOD BOTH TYPES ALWAYS DEFAULT
else !READ DATA FROM EXTERNAL FILES
    open (10, FILE='graphs\v_dss.txt')
    open (17, FILE='graphs\b0_next_dss.txt')
    open (18, FILE='graphs\w_dss.txt')
    do i_b = 1,b_num
        do i_y = 1,y_num
            !WRITE(nout,*) i_b, i_y
            READ(10, '(F15.10, X, F15.10, X, F15.10)') v_matrix(i_b, i_y), v0_matrix(i_b, i_y), v1_matrix(i_y)
            if (v0_matrix(i_b, i_y)< v1_matrix(i_y)) then
                default_decision(i_b, i_y) =2
            else
                default_decision(i_b, i_y) =1
            end if
            READ(18, '(F15.10, X, F15.10, X, F15.10)') w_matrix(i_b, i_y), w0_matrix(i_b, i_y), w1_matrix( i_y)
            READ(17, '(F15.11)') b0_next_matrix(i_b, i_y)
        end do
    end do
end if
CLOSE(10)
CLOSE(17)
CLOSE(18)

call iterate
!call read_data
call simulate
!call compute_initial_guess

call cpu_time(end_time)
WRITE(nout, '(A7, X, A7, X, A7)') 'Hours ', 'Minutes', 'Seconds'
WRITE(nout, '(I7, X, I7, X, I7)') INT((end_time - start_time) / 3600d+0), &
             INT((end_time-start_time)/60d+0 - INT((end_time - start_time) / 3600d+0)*60d+0),&
INT(end_time-start_time - INT((end_time - start_time) / 3600d+0)*3600d+0 - &
INT((end_time-start_time)/60d+0 - INT((end_time - start_time) / 3600d+0)*60d+0)*60d+0)
end program


subroutine read_data
USE param
DOUBLE PRECISION :: b_next
INTEGER :: i_b, i_y, index_b, i

open (101, FILE='graphs\b_grid_dss.txt')
READ(101, '(F12.8)') b_grid(1)

do i=2,b_num
    READ(101, '(F12.8)') b_grid(i)
end do

CLOSE(101)

open (10, FILE='graphs\v_dss.txt')
open (13, FILE='graphs\b_next_dss.txt')

do i_b = 1,b_num
    do i_y = 1,y_num
        READ(10, '(F15.10, X, F15.10, X, F15.10)') v_matrix(i_b, i_y), v0_matrix(i_b, i_y), v1_matrix( i_y)
        if (v0_matrix(i_b, i_y) > v1_matrix( i_y)) then
            default_decision(i_b, i_y) = 1
        else
            default_decision(i_b, i_y) = 2
        end if
        READ(13, '(F15.11)') b_next
        index_b = b_num - 1
        if (ABS(b_next-b_grid(1)) < 1d-8) then
            index_b = 1
        end if
        do i=1,b_num-1
            IF(ABS(b_next-b_grid(i)) < 1d-8) then !(b_next > b_grid(i) .AND. b_next<=b_grid(i+1)) then
                index_b = i
            end if
        end do
        b_next_matrix(i_b, i_y) = index_b
        if (ABS(b_next - b_grid(index_b))>.00001)  then
            WRITE(nout, *)
            pause
        end if
    end do
end do

CLOSE(10)
CLOSE(13)
end subroutine

subroutine simulate
USE param
integer :: period_num, i,j,k, gov_type, random_num, ivalue, sample_num, MAXFN, default_num, default_pol_num,&
           counter, i_y_previous, i_y_current, i_b_current
parameter (sample_num =1, period_num=500001)
DOUBLE PRECISION :: random_matrix(period_num, sample_num, 3), random_vector(1:3*sample_num*period_num), &
                    z(period_num,sample_num), b(period_num+1,sample_num), &
                    q(period_num,sample_num), eps, v_def, v_no_def, b_next, q_fun, qdef, &
                    q_paid, b_zero, c(period_num,sample_num), tb(period_num,sample_num),&
                    y(period_num,sample_num), objective_function
INTEGER :: d(period_num, sample_num), excl(period_num,sample_num),default_potential(period_num, sample_num), liquidity(period_num, sample_num)
EXTERNAL q_fun, objective_function

b_zero = 0d+0  !IF THE COUNTRY IS NOT EXCLUDED AFTER A DEFAULT EPISODE ==> BORROWS AS IF IT STARTS WITH ZERO DEBT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! NOTE !!!!!!
! z DENOTE UNDERLYING SHOCK TO THE ENDOWMENT
! y DENOTE REALIZED ENDOWMENT (EXP(z))
! CODE WAS WRITEN WITH y = SHOCK, SO mean_y ACTUALLY = E(z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



STEP = (b_grid(2) - b_grid(1))*0.5
BOUND = 10*y_sup
XACC = 1d-6
MAXFN = 1000

ivalue=139719
call RNSET(ivalue)       !INITIALIZE THE SEED SO THE SAME RANDOM NUMBERS ARE GENERATED IN EVERY REALIZATION
call DRNUN(2*period_num*sample_num, random_vector)


!First column of random_matrix is used to generate transitory shocks.
!Second column of random_matrix is used to generate type changes.
do j=1,sample_num
    do i=1,period_num
       random_matrix(i,j,1)=random_vector((j-1)*period_num+i)
       random_matrix(i,j,2)=random_vector(period_num*sample_num + (j-1)*period_num+i)
       random_matrix(i,j,3)=random_vector(2*period_num*sample_num + (j-1)*period_num+i)
    end do
end do



open (UNIT=21, FILE="graphs\data_sim_dss.txt", status = 'replace')
open (UNIT=22, FILE="graphs\def_per_dss.txt", status = 'replace')
open (UNIT=23, FILE="graphs\param_dss.txt", status = 'replace')
open (UNIT=24, FILE="graphs\tirar_dss.txt", status = 'replace')
open (UNIT=25, FILE="graphs\def_pol_per_dss.txt", status = 'replace')
open (UNIT=26, FILE="graphs\def_pol_num_dss.txt", status = 'replace')

WRITE(23, '(I10)') period_num-1
WRITE(23, '(I10)') sample_num
CLOSE(23)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
i_default_global = 1  !USED TO INVOKE WHICH CHEBYCHEV MATRIX TO USE
                      !NEED TO MODIFY THIS WHEN DEFAULT AFFECTS OUTPUT REGARDLESS OF THE EXCLUSION STATUS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=1,sample_num   !SOLVE FOR SAMPLE j
    WRITE(nout, *) j
    i_y_previous = INT(y_num)/2
    z(1,j) = y_grid(i_y_previous)
    y(1,j) = EXP(z(1,j))
    b(1,j) = zero
    i_b_current = i_b_zero                   !INDEX OF BOND POSITION IN THE NEXT PERIOD.
    b(2,j) = b_grid(i_b_zero)
    d(1,j) = 1    !NO DEFAULT IN FIRST PERIOD
    excl(1,j) = 1 !NO EXCLUSION IN FIRST PERIOD
    default_num =0d+0  !variable that counts the number of defaults in each sample.
    default_pol_num =0d+0  !variable that counts the number of defaults caused by political change in each sample.
    liquidity(1,j) = 1  !ACCESS TO MARKETS IN THE FIRST PERIOD, NO SS or interest rate shock in the first period
    i_ss_initial = liquidity(1,j)
    do i=2,period_num
    !epsilon = realization of standard gaussian * standard deviation
    !epsilon = realization of standard gaussian * standard deviation
    qdef = 0d+0
    i_y_current = 1
    do WHILE(cdf_matrix(i_y_previous, i_y_current) < random_matrix(i,j,1))
        i_y_current = i_y_current + 1
    end do
    !GROWTH SHOCK
    z(i,j) = y_grid(i_y_current)
    i_y_previous = i_y_current
    i_b_current = i_b_next  !CURRENT INDEX OF ASSET = PREVIOUS INDEX OF SAVING
    y_initial = z(i,j)
    b_initial = b(i,j)
    if (random_matrix(i,j,3) <=trans_matrix(liquidity(i-1,j), 2) ) THEN  !SWITCH FROM STATE liquidity(i-1,j) TO SUDDEN STOP
        liquidity(i,j) = 2
    else
        liquidity(i,j) = 1
    end if
    i_ss_initial = liquidity(i,j)
    if (excl(i-1,j) ==2 ) THEN  !COUNTRY WAS EXCLUDED IN THE PREVIOUS PERIOD
        if (random_matrix(i,j,2) < prob_excl_end) THEN  ! EXCLUSION ENDS THIS PERIOD
            d(i,j) = 1   !COUNTRY DOES NOT DEFAULT
            excl(i,j) = 1
            i_b_next = b_next_matrix(i_b_current, i_y_current)
            b(i+1,j) = b_grid(i_b_next)
            q(i,j) = q_fun(i_b_next, i_y_current)
            y(i,j) = scale_y * EXP(z(i,j))
            c(i,j) = y(i,j) + b(i,j) - q(i,j)*b(i+1,j)
            tb(i,j) = y(i,j) - c(i,j)
        else ! EXCLUSION DOES NOT END THIS PERIOD
            d(i,j) = 1
            excl(i,j) = 2
            i_b_next = i_b_zero
            b(i+1,j) = b_grid(i_b_next)
            q(i,j) = q_fun(i_b_next, i_y_current)
            y(i,j) = scale_y * MIN(lambda, exp(z(i,j)))
            c(i,j) = y(i,j) + b(i,j) - q(i,j)*b(i+1,j)
            tb(i,j) = y(i,j) - c(i,j)
        end if
    else !COUNTRY WAS NOT EXCLUDED IN THE PREVIOUS PERIOD
        !WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8)') b_initial, v_no_def, v_def, 100*random_matrix(i,j,1)
        if (default_decision(i_b_current, i_y_current) > 1) then  !COUNTRY DEFAULTS
            d(i,j) = 2
            excl(i,j) = 2
            i_b_next = i_b_zero
            b(i+1,j) = b_grid(i_b_next)
            q(i,j) = q_fun(i_b_next, i_y_current)  !q_paid
            y(i,j) = scale_y * MIN(lambda, exp(z(i,j)))
            c(i,j) = y(i,j)
            tb(i,j) = y(i,j) - c(i,j)
            qdef = q_fun(b0_next_matrix_indice(i_b_current, i_y_current),i_y_current)
            WRITE(22, '(I7)') i-1 !NEED TO SUBSTRACT 1. REASON: files start saving data on period 2
            default_num = default_num + 1
        else
            d(i,j) = 1   !COUNTRY DOES NOT DEFAULT
            excl(i,j) = 1
            i_b_next = b_next_matrix(i_b_current, i_y_current)
            b(i+1,j) = b_grid(i_b_next)
            !Compute the bond price paid when an amount b_next is issued.
            !the bond price depends on the current default decision and output (both will affect output tomorrow)
            !It also depends on the current type in power.

            q(i,j) = q_fun(i_b_next, i_y_current)  !q_paid
            y(i,j) = scale_y * EXP(z(i,j))
            c(i,j) = y(i,j) + b(i,j) - q(i,j)*b(i+1,j)
            tb(i,j) = y(i,j) - c(i,j)
        end if
    end if
    WRITE(21, '(F12.8, X, F12.8, X, F12.8, X, F12.8,  X, F12.8, X, F12.8, X, I3, X, I3)') &
    LOG(y(i,j)), b(i,j), q(i,j), qdef, LOG(c(i,j)), tb(i,j)/y(i,j), excl(i,j), d(i,j)
    end do
    WRITE(26, '(I6, X, I6)') default_num, default_pol_num
end do
CLOSE(21)
CLOSE(22)
CLOSE(24)
CLOSE(25)
CLOSE(26)
end subroutine