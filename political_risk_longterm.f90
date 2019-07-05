!Last change:  May 3, 2016   10:53, YKO
module param
DOUBLE PRECISION :: beta, sigma, r, b_inf, b_sup, y_inf, y_sup, eps_inf, eps_sup, std_eps,&
                    rho, mean_y, std_y, pi_number, lambda, zero, width, prob_excl_end, delta, scale_y, &
                    cdf_inf, cdf_sup, coupon, alpha0, alpha1, d0, d1, prob_excl_end1, phi, prob_office, theta, share

parameter (beta = 0.96d+0, sigma = 2d+0, r = 0.01d+0, &
           std_eps = 0.025d+0, rho = 0.945d+0, mean_y = -0.5d+0*std_eps**2,  pi_number = 3.1415926535897932d+0, prob_office = 0.85d+0, &
          lambda = 0.969d+0, zero = 0d+0,  width = 1.5d+0, prob_excl_end = 0.282d+0,  delta = 0.0375d+0, & !0.045d+0, &
          cdf_inf = 3.167d-5, cdf_sup = 1d+0 - 3.167d-5, scale_y = 1d+0, &
          d0 = -0.69d+0, d1 = 1.01d+0, prob_excl_end1 = 0.0d+0, phi = 0.0d+0, theta = 0.62d+0, share = (theta/(1d+0 - theta))**(1d+0/sigma))

!y = log(income)
!sigma = coefficient of relative risk aversion
!r = international interest rate
!rho = coefficient of autocorrelation in income
!std_eps = standard deviation of innovation to income shocks
!mean_y = unconditional mean of the growth trend shock
!std_y = standard deviation of growth trend shock


INTEGER :: b_num, y_num, eps_num, quad_num_v, quad_num_q, nout, i_y_global, b_num_finer, y_num_finer, &
           i_b_global, i_default_global,  cdf_num
parameter (b_num= 25, y_num = 25, quad_num_v = 50, eps_num = quad_num_v, quad_num_q = 50, cdf_num = 50, b_num_finer = 50, y_num_finer = 50)

DOUBLE PRECISION :: b_grid(1:b_num), y_grid(1:y_num), default_grid(1:2), indicator_tirar,indicator_tirar1,&
                    y_initial, b_initial, b_global, eps_grid(1:eps_num), cdf_grid(cdf_num), BREAK_eps(cdf_num), &
                    CSCOEF_eps(4,cdf_num), quad_w_v(1:quad_num_v), quad_x_v(1:quad_num_v), quad_w_hermite(1:eps_num),&
                    quad_x_hermite(1:eps_num), quad_w_q(1:quad_num_q), quad_x_q(1:quad_num_q), counter, y_grid_finer(1:y_num_finer), b_grid_finer(1:b_num_finer),&
                    y_global

DOUBLE PRECISION, DIMENSION(b_num, y_num) :: v_matrix, w_matrix, default_decision, b_next_matrix, v_matrix_nodil
DOUBLE PRECISION, DIMENSION(b_num, y_num) :: v0_matrix, v1_matrix, b0_next_matrix, b1_next_matrix, q_matrix, q_matrix_nodef,&
                                             q_matrix_nodef_rn, q_matrix_rn, v_excl_matrix, b_adj_matrix, w0_matrix, w1_matrix, w_excl_matrix

DOUBLE PRECISION, DIMENSION(y_num, b_num) :: break_matrix, break_matrix_q, break_matrix_v1, break_matrix_w, break_matrix_w1, break_matrix_w0
DOUBLE PRECISION, DIMENSION(y_num, 4, b_num) :: coeff_matrix, coeff_matrix_q, coeff_matrix_v1, coeff_matrix_w, coeff_matrix_w1, coeff_matrix_w0

DOUBLE PRECISION, DIMENSION (b_num_finer, y_num_finer) :: EV_matrix, q_menu_matrix, EW_matrix
DOUBLE PRECISION, DIMENSION (b_num_finer, y_num_finer) :: EV_excl_matrix, EW_excl_matrix

DOUBLE PRECISION, DIMENSION(b_num_finer, y_num_finer) :: break_matrix_EV, break_matrix_q_menu, break_matrix_EV_excl, break_matrix_EW, break_matrix_EW_excl
DOUBLE PRECISION, DIMENSION(4, b_num_finer, y_num_finer) :: coeff_matrix_EV, coeff_matrix_q_menu, coeff_matrix_EV_excl, coeff_matrix_EW, coeff_matrix_EW_excl


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
!WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)')  c, cons1, cons2, u_funw
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
DOUBLE PRECISION :: dos, DNORDF, delta_y, prob_mass, b_inf0, b_sup0, b_inf1, b_sup1
INTEGER :: b_num_half
EXTERNAL :: DNORDF

std_y = std_eps/ SQRT(1 - rho**2)

coupon = (r+delta)/(1d+0 + r)
y_inf = mean_y - 4*std_y
y_sup = mean_y + 4*std_y
delta_y = 4*std_y / (y_num - 1d+0)

b_inf = -3.0d+0 !* EXP(y_inf)*lambda / (1-recovery)
b_sup =  -0.0001d+0

open (11, FILE='graphs\y_grid.txt',STATUS='replace')


do i=1,b_num
   b_grid(i) = b_inf + (b_sup - b_inf) * (i-1) / (b_num - 1)
   !WRITE(10, '(F12.8)') b_grid(i)
end do

do i=1,y_num
   y_grid(i) = y_inf + (y_sup - y_inf) * (i-1) / (y_num - 1)
   WRITE(11, '(F12.8)') y_grid(i)
end do

do i=1,cdf_num
   cdf_grid(i) = cdf_inf + (cdf_sup - cdf_inf) * (i-1) / (cdf_num - 1)
end do

CLOSE(11)


default_grid(1) = 0.0d+0
default_grid(2) = 1.0d+0

eps_inf = -4*std_eps
eps_sup = 4*std_eps

do i=1,eps_num
   !GRID WHEN USING GAUSS-HERMITE QUADRATURE POINTS TO COMPUTE EXPECTED VALUES
   eps_grid(i) = SQRT(dos) * std_eps* quad_x_hermite(i)

   !GRID WHEN USING GAUSS-LEBESGUE QUADRATURE POINTS TO COMPUTE EXPECTED VALUES
   !NEED TO DIFFERENTIATE BECAUSE THE FIRST GRID ASSIGNS TO MUCH WEIGHT ON POINTS
   !WITH 'ZERO' DENSITY. THIS IS INEFFICIENT IF THE INTEGRALS ARE COMPUTED USING
   !GAUSS-LEBESGUE QUADRATURE POINTS.
!   eps_grid(i) = 4* std_eps* quad_x(i)
end do
do i=1,b_num_finer
   b_grid_finer(i) = b_inf + (b_sup - b_inf)*(i-1)/(b_num_finer - 1)
end do
do i=1,y_num_finer
   y_grid_finer(i) = y_inf + (y_sup - y_inf) * (i-1) / (y_num_finer - 1)
end do

end subroutine



!COMPUTE QUADRATURE POINTS AND WEIGHTS USING LEGENDRE AND GAUSSIAN QUADRATURE RULES
!THEY ARE USED TO COMPUTE NUMERICAL INTEGRALS
subroutine quadrature
USE param
INTEGER :: N_v, N_q, N, IWEIGH, IWEIGH4, NFIX
parameter(N_v = quad_num_v, N_q = quad_num_q, N = eps_num)
DOUBLE PRECISION:: ALPHA, BETA1, QW_v(1:N_v), QX_v(1:N_v), QW_q(1:N_q), QX_q(1:N_q),QW(1:N), QX(1:N), QXFIX(2)
PARAMETER(ALPHA=0, BETA1=0, IWEIGH=1, IWEIGH4=4)
external DGQRUL

!quad_w = weights using Gauss legendre quadrature rule
!quad_x = points using Gauss legendre quadrature rule
!quad_w_hermite = weights using Gauss Hermite quadrature rule
!quad_x_hermite = points using Gauss Hermite quadrature rule

NFIX = 0
CALL DGQRUL(N_v,IWEIGH, ALPHA, BETA1, NFIX, QXFIX, QX_v, QW_v)

quad_w_v=QW_v
quad_x_v=QX_v

CALL DGQRUL(N_q,IWEIGH, ALPHA, BETA1, NFIX, QXFIX, QX_q, QW_q)

quad_w_q=QW_q
quad_x_q=QX_q


CALL DGQRUL(N,IWEIGH4, ALPHA, BETA1, NFIX, QXFIX, QX, QW)
quad_x_hermite = QX
quad_w_hermite = QW


end subroutine

!COMPUTES THE DIFFERENCE BETWEEN THE VALUE FUNCTION UNDER NO DEFAULT AND DEFAULT
!NEED b_global TO BE DEFINED BEFORE.
DOUBLE PRECISION function dif_fun(y)
USE param
DOUBLE PRECISION, INTENT(IN) :: y
DOUBLE PRECISION :: v0_fun, v1_fun
EXTERNAL v0_fun, v1_fun


!dif_fun = interpolate(b_global, y, v0_matrix) - interpolate(b_global, y, v1_matrix)
dif_fun = v0_fun(b_global, y) - v1_fun(b_global, y)
!write(nout, '(F12.8, X, F12.8, X, F12.8)') b_global, v0_fun(b_global, y), v1_fun(b_global,y)
end

!FUNCTION USED TO COMPUTE THE MINIMUM INCOME FOR WHICH THE CURRENT PARTY IN POWER DOES NOT DEFAULT
!b = outstanding debt
DOUBLE PRECISION function y_fun(b)
USE param
DOUBLE PRECISION, INTENT(IN) :: b
INTEGER :: MAXFN, num_tirar
PARAMETER(num_tirar = 100)
DOUBLE PRECISION :: dif_fun, ERRABS, ERRREL, left, right, y_max, y_min, dif_right, dif_left,&
                    y_tirar(num_tirar), tirar
EXTERNAL dif_fun, DZBREN

MAXFN=1000
ERRREL = 1D-10
ERRABS = 1D-10


y_max = rho*y_initial + (1-rho)*mean_y + width*eps_sup
y_min = rho*y_initial + (1-rho)*mean_y + width*eps_inf

b_global = b     !b_global IS USED IN dif_fun TO FIX THE VALUE OF b (OUTSTANDING DEBT)

!open (120, FILE='graphs\tirar1.txt',STATUS='replace')
!do i=1,num_tirar
!    tirar = y_min + (y_max - y_min) * (i-1d+0)/(num_tirar-1)
!    left = dif_fun(tirar)
!end do
!CLOSE(120)


dif_right = dif_fun(y_max)
dif_left  = dif_fun(y_min)

!1) DETERMINE WHETHER THERE IS AN INTERIOR ROOT OR NOT
!   a) IF THERE IS AN INTERIOR ROOT, SPECIFY WHETHER THE difference function IS INCREASING IN g OR NOT.
!   b) IF THERE IS NO INTERIOR ROOT, DETERMINE WHETHER THE GOV'T ALWAYS OR NEVER DEFAULTS


if (dif_right*dif_left<0) then  !THERE IS AN INTERIOR ROOT
   left =  y_min
   right = y_max
   CALL DZBREN (dif_fun, ERRABS, ERRREL, left, right, MAXFN)
   y_fun = right
else if (dif_left>0) then
      y_fun = y_min
else !dif_right<0
      y_fun = y_max
end if
end


DOUBLE PRECISION function q_fun(b_next, y)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_next, y
DOUBLE PRECISION :: prob_no_def_tomorrow, y_current_type, DNORDF, y_fun, standarized_inf, standarized_sup,&
                    y_threshold, y_min, y_max, exp_q, interpolate, scalar, y_next, q_paid_fun, acum_int, &
                    DCSVAL, cdf, cdf_threshold, DNORIN, m_next, exp_m, exp_q_rn, q_paid_fun_rn, prop, gamma
EXTERNAL DNORDF, y_fun, interpolate, q_paid_fun, DCSVAL, DNORIN, q_paid_fun_rn


!y_current_type = STANDARIZED INCOME THRESHOLDS IF THE TYPE DOES NOT CHANGE TOMORROW

exp_m = 0d+0
exp_q = 0d+0
scalar = 1d+0 / SQRT(pi_number)
NINTV = cdf_num - 1

if (b_next>0) then
    !NO SAVING
    q_fun = 1000
ELSEIF(b_next >=b_inf) THEN !SET POSITIVE PRICES ONLY FOR b' THAT ARE ABOVE THE MINIMUM VALUE.
                            !(AVOID OPTIMAL VALUES OUTSIDE THE RANGE)
    y_threshold = y_fun(b_next)
    y_current_type = (y_threshold - rho * y - (1-rho)*mean_y) / std_eps
    cdf_threshold = DNORDF(y_current_type)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EN SIGUIENTES LINEAS INICIALIZO EL exp_m A CERO, DONDE exp_m ES EL VALOR ESPERADO DE m_next, Y QUE SALE
! DE INTEGRAR m_next IGUAL A COMO SE HACE CON exp_q .	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (cdf_threshold < cdf_inf) then !THRESHOLD TOMORROW IS LOW
        prob_no_def_tomorrow = 1d+0
        !WRITE(nout, *) 'prob no def tomorrow = ', prob_no_def_tomorrow
        do i=1,quad_num_q
            cdf = 0.5*(quad_x_q(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
            y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
            m_next = EXP(-r)
            exp_m = exp_m+ m_next * coupon * 0.5* (cdf_sup - cdf_inf)*quad_w_q(i)
            exp_q = exp_q + (1d+0 - delta) * m_next * q_paid_fun(b_next, y_next) * 0.5* (cdf_sup - cdf_inf)*quad_w_q(i)
        end do
            !NEED TO ADJUST EXPECTATION TO TAKE INTO ACCOUNT THAT THE
            !ACTUAL PROBABILITY MASS BETWEEN y_threshold and infty IS UNDERESTIMATED USING LEBESGE-..
            exp_q = exp_q / (cdf_sup - cdf_inf)  !NEED TO ADJUST EXPECTATION TO TAKE INTO ACCOUNT THAT THE
            exp_m = exp_m / (cdf_sup - cdf_inf)  !!NUEVO
    elseif (cdf_threshold > cdf_sup) then !THRESHOLD TOMORROW IS HIGH
        prob_no_def_tomorrow = 0d+0
        exp_q = 0d+0
        exp_m = 0d+0
    else
        prob_no_def_tomorrow = (1d+0 - cdf_threshold)
        !WRITE(nout, *) 'prob no def tomorrow = ', prob_no_def_tomorrow
        do i=1,quad_num_q
            !compute integral for y > y_threshold
            cdf = 0.5*(quad_x_q(i) + 1)*(cdf_sup - cdf_threshold) + cdf_threshold
            y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
            m_next = EXP(-r )
            exp_m = exp_m+  m_next * coupon * 0.5* (cdf_sup - cdf_threshold)*quad_w_q(i)
            exp_q = exp_q + (1d+0 - delta) * m_next * q_paid_fun(b_next, y_next) * 0.5* (cdf_sup - cdf_threshold)*quad_w_q(i)
        end do
        !NEED TO ADJUST EXPECTATION TO TAKE INTO ACCOUNT THAT
        !THE ACTUAL PROBABILITY MASS BETWEEN y_threshold and infty IS UNDERESTIMATED USING LEBESGE-..
        exp_q = exp_q / (cdf_sup - cdf_inf)  !NEED TO ADJUST EXPECTATION TO TAKE INTO ACCOUNT THAT
        exp_m = exp_m / (cdf_sup - cdf_inf)  !NEED TO ADJUST EXPECTATION TO TAKE INTO ACCOUNT THAT
    end if
q_fun = MIN(MAX(exp_m + exp_q, 0d+0),1d+0/(1d+0+r)) 
end if
!WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') b_next, y, q_fun, exp_m, exp_q
!pause
end



    
DOUBLE PRECISION function objective_excl(b)
USE param
DOUBLE PRECISION, INTENT(IN) :: b
INTEGER :: other_type, i,j,h, t, d, NINTV
DOUBLE PRECISION :: u_fun, q_fun, acum, exp_v_next, value_next, y_next, y_fun, scalar, DCSVAL, &
                    y_threshold, y_max, y_min, y_next1, acum2,acum1, g, q, DNORDF, v0_fun, v1_fun, &
                    value_next_tirar, acum_int, acum_int_lo, acum_int_hi, borrowing, interpolate, cdf, cdf1,&
                    cdf_threshold, DNORIN, q_risk_free, output, v_excl_fun, b_next, EV_excl_fun, w0_fun, w1_fun, exp_w_next

EXTERNAL u_fun, q_fun, DNORDF, y_fun, interpolate, v0_fun, v1_fun, DCSVAL, DNORIN, v_excl_fun, EV_excl_fun, w0_fun, w1_fun

d  = i_default_global
b_next  = MAX(MIN(b, b_sup), b_inf)
exp_v_next = EV_excl_fun(b_next, y_initial)
!exp_w_next = prob_excl_end*w0_fun(b_next, y_initial) + (1-prob_excl_end)*w1_fun(b_next, y_initial)

!output = scale_y * MIN(lambda, exp(y_initial))
output = (1d+0 - default_grid(d)) * EXP(y_initial) + default_grid(d) * (EXP(y_initial) *(1d+0 - d0) - d1*EXP(y_initial)**2d+0)
!output = MIN(lambda, exp(y_initial))
objective_excl = - u_fun(output) - beta *exp_v_next !(prob_office*exp_v_next + (1-prob_office)*exp_w_next)
end
    
DOUBLE PRECISION function objective_excl_w(b)
USE param
DOUBLE PRECISION, INTENT(IN) :: b
INTEGER :: other_type, i,j,h, t, d, NINTV
DOUBLE PRECISION :: u_fun_w, q_fun, acum, exp_v_next, value_next, y_next, y_fun, scalar, DCSVAL, &
                    y_threshold, y_max, y_min, y_next1, acum2,acum1, g, q, DNORDF, v0_fun, v1_fun, w0_fun, w1_fun, &
                    value_next_tirar, acum_int, acum_int_lo, acum_int_hi, borrowing, interpolate, cdf, cdf1,&
                    cdf_threshold, DNORIN, q_risk_free, output, v_excl_fun, b_next, EV_excl_fun, exp_w_next, EW_excl_fun

EXTERNAL u_fun_w, q_fun, DNORDF, y_fun, interpolate, v0_fun, v1_fun, DCSVAL, DNORIN, v_excl_fun, EV_excl_fun, EW_excl_fun, w0_fun, w1_fun

d  = i_default_global
b_next  = MAX(MIN(b, b_sup), b_inf)
exp_v_next = EW_excl_fun(b_next, y_initial)
!exp_w_next = prob_excl_end*w0_fun(b_next, y_initial) + (1-prob_excl_end)*w1_fun(b_next, y_initial)

output = (1d+0 - default_grid(d)) * EXP(y_initial) + default_grid(d) * (EXP(y_initial) *(1d+0 - d0) - d1*EXP(y_initial)**2d+0)

!output = scale_y * MIN(lambda, exp(y_initial))

!output = MIN(lambda, exp(y_initial))
objective_excl_w = - u_fun_w(output) - beta*exp_v_next !*(prob_office*exp_w_next + (1-prob_office)*exp_v_next)
!WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') b_next, exp_v_next, w0_fun(b_next, y_initial), w1_fun(b_next, y_initial), objective_excl_w
end    

!Compute the objective function in the Belman equation
DOUBLE PRECISION function objective_function(b)
USE param
DOUBLE PRECISION, INTENT(IN) :: b
INTEGER :: other_type, i,j,h, t, d, NINTV
DOUBLE PRECISION :: u_fun, q_fun, acum, exp_v_next, value_next, y_next, y_fun, scalar, DCSVAL, &
                    y_threshold, y_max, y_min, y_next1, acum2,acum1, g, q, DNORDF, v0_fun, v1_fun, &
                    value_next_tirar, acum_int, acum_int_lo, acum_int_hi, borrowing, interpolate, cdf, cdf1,&
                    cdf_threshold, DNORIN, b_nodil, output, bnext, EV_fun, q_menu_fun, gamma, exp_w_next, EW_fun

EXTERNAL u_fun, q_fun, DNORDF, y_fun, interpolate, v0_fun, v1_fun, DCSVAL, DNORIN, EV_fun, q_menu_fun, EW_fun

bnext = max(min(b, b_sup), b_inf)
exp_v_next = EV_fun(bnext, y_initial)
!exp_w_next = EW_fun(bnext, y_initial)
q = q_menu_fun(bnext, y_initial)
d  = i_default_global

borrowing =  bnext - b_initial*(1d+0-default_grid(d))*(1d+0 - delta)
!output = (1d+0 - default_grid(d)) * EXP(y_initial) + default_grid(d) * MIN(lambda, exp(y_initial))
!output = scale_y * exp(y_initial)
output = (1d+0 - default_grid(d)) * EXP(y_initial) + default_grid(d) * (EXP(y_initial) *(1d+0 - d0) - d1*EXP(y_initial)**2d+0)
objective_function = - u_fun(output + coupon*b_initial*(1d+0 -default_grid(d)) - borrowing*q) - beta *exp_v_next !(prob_office*exp_v_next + (1-prob_office)*exp_w_next)

end

!Compute the objective function in the Belman equation
DOUBLE PRECISION function objective_function_w(b)
USE param
DOUBLE PRECISION, INTENT(IN) :: b
INTEGER :: other_type, i,j,h, t, d, NINTV
DOUBLE PRECISION :: u_fun_w, q_fun, acum, exp_v_next, value_next, y_next, y_fun, scalar, DCSVAL, &
                    y_threshold, y_max, y_min, y_next1, acum2,acum1, g, q, DNORDF, v0_fun, v1_fun, &
                    value_next_tirar, acum_int, acum_int_lo, acum_int_hi, borrowing, interpolate, cdf, cdf1,&
                    cdf_threshold, DNORIN, b_nodil, output, bnext, EV_fun, q_menu_fun, gamma, exp_w_next, EW_fun

EXTERNAL u_fun_w, q_fun, DNORDF, y_fun, interpolate, v0_fun, v1_fun, DCSVAL, DNORIN, EV_fun, q_menu_fun, EW_fun

bnext = max(min(b, b_sup), b_inf)
!exp_v_next = EV_fun(bnext, y_initial)
exp_v_next = EW_fun(bnext, y_initial)
q = q_menu_fun(bnext, y_initial)
d  = i_default_global

borrowing =  bnext - b_initial*(1d+0-default_grid(d))*(1d+0 - delta)
!output = (1d+0 - default_grid(d)) * EXP(y_initial) + default_grid(d) * MIN(lambda, exp(y_initial))
!output = (1d+0 - default_grid(d)) * EXP(y_initial) + default_grid(d) * (EXP(y_initial) *(1d+0 - d0) - d1*EXP(y_initial)**2d+0)
output = exp(y_initial)
objective_function_w = - u_fun_w(output + coupon*b_initial*(1d+0 -default_grid(d)) - borrowing*q) - beta *exp_v_next!(prob_office*exp_w_next + (1-prob_office)*exp_v_next)

end

subroutine optimize(b_next, v_value)
USE param
integer :: MAXFN, t, d, index_opt, index_vector(1), num, num_tirar
parameter (num=25, num_tirar = 500)
DOUBLE PRECISION :: b_next, v_value, objective_function, new_value, q_fun, old_value, b, u_fun,q,vector(b_num),&
                    b_next_grid(num), STEP, BOUND, XACC, b_next_guess, b_tirar, b_next_inf, b_next_sup, b_tirar_sup, b_tirar_inf
EXTERNAL objective_function, q_fun, u_fun, DUVMIF

!if (indicator_tirar>0) then
!   open (117, FILE='graphs\tirar1.txt',STATUS='replace')
!end if

b_next_inf = b_inf
b_next_sup = b_sup

old_value = 10d+3
index_opt=1
do i=1,num

   b_next_grid(i) = b_next_inf + (b_next_sup - b_next_inf) * (i-1)/ (num - 1)
   new_value = objective_function(b_next_grid(i))
!   WRITE(nout, *) new_value
   if (new_value <= old_value) then
      index_opt = i
      b_next_guess = b_next_grid(i)
      old_value = new_value
  end if
end do
STEP = (b_next_grid(2) - b_next_grid(1))*0.5
BOUND = 10*EXP(y_sup)
XACC = 1d-6
MAXFN = 1000


b_tirar_sup = b_next_sup - 1.0d-6
b_tirar_inf = b_next_inf + 1.0d-6


if (index_opt == 1 .AND. objective_function(b_tirar_inf) >= objective_function(b_next_grid(index_opt))) then
   b_next = b_next_guess
   v_value = -old_value

elseif(index_opt == 1 .AND. objective_function(b_tirar_inf) < objective_function(b_next_grid(index_opt))) then
       !THE MAXIMUM MAY BE CLOSE TO MINIMUM GRID POINT BUT NOT AT MINIMUM GRID POINT
       b_next_guess = b_tirar_inf
       call DUVMIF(objective_function, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = -objective_function(b_next)
elseif (index_opt == num .AND. objective_function(b_tirar_sup) >= objective_function(b_next_grid(index_opt))) then
       b_next = b_next_guess
       v_value = -old_value
elseif(index_opt == num .AND. objective_function(b_tirar_sup) < objective_function(b_next_grid(index_opt))) then
       !THE MAXIMUM MAY BE CLOSE TO MAXIMUM GRID POINT BUT NOT AT MAXIMUM GRID POINT
       b_next_guess = b_tirar_sup
       call DUVMIF(objective_function, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = -objective_function(b_next)
else
   if (ABS(old_value - objective_function(b_next_grid(index_opt-1)))/(b_next_grid(2)-b_next_grid(1)) < 1d-10 .OR.  &
       ABS(old_value - objective_function(b_next_grid(index_opt+1)))/(b_next_grid(2)-b_next_grid(1)) < 1d-10) then
       b_next = b_next_guess
       v_value = -old_value
   else
      call DUVMIF(objective_function, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = -objective_function(b_next)
       indicator_tirar1 =0
   end if
end if

10 end subroutine


DOUBLE PRECISION function interpolate(b, y, matrix)
USE param
DOUBLE PRECISION, INTENT(IN) :: b, y, matrix(b_num, y_num)
DOUBLE PRECISION ::  slope, weight, acum, ratio_b, ratio_y
INTEGER :: index_b, index_y, i_b, i_y

index_y = (MAX(MIN(INT((y_num-1)*(y-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))
index_b = (MAX(MIN(INT((b_num-1)*(b-b_grid(1))/(b_grid(b_num)-b_grid(1)))+1,b_num-1), 1))


ratio_b = (b - b_grid(index_b)) / (b_grid(index_b+1) - b_grid(index_b))
!ratio_y = MIN(MAX((y - y_grid(index_y)) / (y_grid(index_y+1) - y_grid(index_y)),0d+0), 1d+0)
ratio_y = (y - y_grid(index_y)) / (y_grid(index_y+1) - y_grid(index_y))
acum=0

do i_b=0,1
   do i_y=0,1
        weight = ((1-i_b)*(1 - ratio_b) + i_b*ratio_b )* ((1-i_y)*(1 - ratio_y) + i_y*ratio_y )
        acum = acum + matrix(index_b + i_b, index_y + i_y)*weight
   end do
end do


interpolate = acum
    end


DOUBLE PRECISION function interpolate_finer(b, y, matrix)
USE param
DOUBLE PRECISION, INTENT(IN) :: b, y, matrix(b_num_finer, y_num_finer)
DOUBLE PRECISION ::  slope, weight, acum, ratio_b, ratio_y
INTEGER :: index_b, index_y, i_b, i_y

index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))
index_b = (MAX(MIN(INT((b_num_finer-1)*(b-b_grid_finer(1))/(b_grid_finer(b_num_finer)-b_grid_finer(1)))+1,b_num_finer-1), 1))


ratio_b = (b - b_grid_finer(index_b)) / (b_grid_finer(index_b+1) - b_grid_finer(index_b))
!ratio_y = MIN(MAX((y - y_grid(index_y)) / (y_grid(index_y+1) - y_grid(index_y)),0d+0), 1d+0)
ratio_y = (y - y_grid_finer(index_y)) / (y_grid_finer(index_y+1) - y_grid_finer(index_y))
acum=0

do i_b=0,1
   do i_y=0,1
        weight = ((1-i_b)*(1 - ratio_b) + i_b*ratio_b )* ((1-i_y)*(1 - ratio_y) + i_y*ratio_y )
        acum = acum + matrix(index_b + i_b, index_y + i_y)*weight
   end do
end do


interpolate_finer = acum
end
    
DOUBLE PRECISION function v0_fun(b, y)
USE param
DOUBLE PRECISION, INTENT(IN) :: b, y
DOUBLE PRECISION ::  slope, DCSVAL, value_left, value_right, interpolate
INTEGER :: index_b, NINTV
EXTERNAL DCSVAL, interpolate

NINTV = b_num - 1

! THIS IS FOR SPLINE
index_y = (MAX(MIN(INT((y_num-1)*(y-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))
value_left  = DCSVAL(b, NINTV, break_matrix(index_y,:), coeff_matrix(index_y,:,:))
value_right = DCSVAL(b, NINTV, break_matrix(index_y+1,:), coeff_matrix(index_y+1,:,:))
slope = (value_right - value_left)/ (y_grid(index_y+1) - y_grid(index_y))
v0_fun = value_left + slope * (y - y_grid(index_y))

! INSTEAD, WE COULD HAVE JUST LINEAR INTERPOLATION
!v0_fun = interpolate(b,y,v0_matrix)
end

DOUBLE PRECISION function w0_fun(b, y)
USE param
DOUBLE PRECISION, INTENT(IN) :: b, y
DOUBLE PRECISION ::  slope, DCSVAL, value_left, value_right, interpolate
INTEGER :: index_b, NINTV
EXTERNAL DCSVAL, interpolate

NINTV = b_num - 1

! THIS IS FOR SPLINE
index_y = (MAX(MIN(INT((y_num-1)*(y-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))
value_left  = DCSVAL(b, NINTV, break_matrix_w0(index_y,:), coeff_matrix_w0(index_y,:,:))
value_right = DCSVAL(b, NINTV, break_matrix_w0(index_y+1,:), coeff_matrix_w0(index_y+1,:,:))
slope = (value_right - value_left)/ (y_grid(index_y+1) - y_grid(index_y))
w0_fun = value_left + slope * (y - y_grid(index_y))

! INSTEAD, WE COULD HAVE JUST LINEAR INTERPOLATION
!w0_fun = interpolate(b,y,w0_matrix)
end    
    
DOUBLE PRECISION function q_paid_fun(b, y)
USE param
DOUBLE PRECISION, INTENT(IN) :: b, y
DOUBLE PRECISION ::  slope, DCSVAL, value_left, value_right, interpolate, b_local
INTEGER :: index_b, NINTV
EXTERNAL DCSVAL, interpolate

b_local = min(max(b, b_inf+1d-12), b_sup-1d-12)
NINTV = b_num - 1
index_y = (MAX(MIN(INT((y_num-1)*(y-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))
value_left  = DCSVAL(b_local, NINTV, break_matrix_q(index_y,:), coeff_matrix_q(index_y,:,:))
value_right = DCSVAL(b_local, NINTV, break_matrix_q(index_y+1,:), coeff_matrix_q(index_y+1,:,:))
slope = (value_right - value_left)/ (y_grid(index_y+1) - y_grid(index_y))

if (slope<0) then 
    q_paid_fun = value_left
else
    q_paid_fun = min(value_left + slope * (y - y_grid(index_y)), EXP(-r))
end if
!WRITE(nout, '(F12.8, X, F12.8, X, F12.8)') value_left_left, value_left, value_right
!pause
!q_paid_fun = interpolate (b_local,y,q_matrix_nodef)
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION function v1_fun(b, y)
USE param
DOUBLE PRECISION, INTENT(IN) ::b, y
DOUBLE PRECISION ::  slope, value_left, value_right, DCSVAL, interpolate
INTEGER :: index_y
EXTERNAL DCSVAL, interpolate

index_y = (MAX(MIN(INT((y_num-1)*(y-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))

value_left  = v1_matrix(1,index_y)
value_right = v1_matrix(1,index_y+1)
!value_left  = DCSVAL(b, NINTV, break_matrix_v1(index_y,:), coeff_matrix_v1(index_y,:,:))
!value_right = DCSVAL(b, NINTV, break_matrix_v1(index_y+1,:), coeff_matrix_v1(index_y+1,:,:))
slope = (value_right - value_left)/ (y_grid(index_y+1) - y_grid(index_y))
v1_fun = value_left + slope * (y - y_grid(index_y))

!value_left = interpolate(b,y,coeff_matrix(index_y,:,:))
!value_right = interpolate(b,y,coeff_matrix(index_y+1,:,:)) ! THIS SHOULD DO THE LINEAR INTERPOLATION
!slope = (value_right - value_left)/ (y_grid(index_y+1) - y_grid(index_y))
!v1_fun = value_left + slope * (y - y_grid(index_y))
end
    
DOUBLE PRECISION function w1_fun(b, y)
USE param
DOUBLE PRECISION, INTENT(IN) ::b, y
DOUBLE PRECISION ::  slope, value_left, value_right, DCSVAL, interpolate
INTEGER :: index_y
EXTERNAL DCSVAL, interpolate

index_y = (MAX(MIN(INT((y_num-1)*(y-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))

value_left  = w1_matrix(1,index_y)
value_right = w1_matrix(1,index_y+1)
!value_left  = DCSVAL(b, NINTV, break_matrix_v1(index_y,:), coeff_matrix_v1(index_y,:,:))
!value_right = DCSVAL(b, NINTV, break_matrix_v1(index_y+1,:), coeff_matrix_v1(index_y+1,:,:))
slope = (value_right - value_left)/ (y_grid(index_y+1) - y_grid(index_y))
w1_fun = value_left + slope * (y - y_grid(index_y))

!value_left = interpolate(b,y,coeff_matrix(index_y,:,:))
!value_right = interpolate(b,y,coeff_matrix(index_y+1,:,:)) ! THIS SHOULD DO THE LINEAR INTERPOLATION
!slope = (value_right - value_left)/ (y_grid(index_y+1) - y_grid(index_y))
!w1_fun = value_left + slope * (y - y_grid(index_y))
end    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NUEVO!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION function v_excl_fun(y)
USE param
DOUBLE PRECISION, INTENT(IN) :: y
DOUBLE PRECISION ::  slope, value_left, value_right
INTEGER :: index_y


index_y = (MAX(MIN(INT((y_num-1)*(y-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))

value_left  = v_excl_matrix(1,index_y)
value_right = v_excl_matrix(1,index_y+1)

slope = (value_right - value_left)/ (y_grid(index_y+1) - y_grid(index_y))

v_excl_fun = value_left + slope * (y - y_grid(index_y))
end

DOUBLE PRECISION function q_menu_fun(b, y)
USE param
DOUBLE PRECISION, INTENT(IN) :: b, y
DOUBLE PRECISION ::  slope, DCSVAL, value_left, value_right, interpolate_finer, b_local
INTEGER :: index_y, NINTV
EXTERNAL DCSVAL, interpolate_finer

!!RESTRICT TO b POSITION INSIDE THE GRID SPACE. OTHERWISE, INTERPOLATION PROCEDURE FAILS
b_local = min(max(b, b_inf+1d-12), b_sup-1d-12)

NINTV = b_num_finer - 1
index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))
!!
!!value_left_left = DCSVAL(b, NINTV, break_matrix_q(index_y-1,:), coeff_matrix_q(index_y-1,:,:))
value_left  = DCSVAL(b_local, NINTV, break_matrix_q_menu(:, index_y), coeff_matrix_q_menu(:,:, index_y))
value_right = DCSVAL(b_local, NINTV, break_matrix_q_menu(:, index_y+1), coeff_matrix_q_menu(:,:, index_y+1))
slope = (value_right - value_left)/ (y_grid_finer(index_y+1) - y_grid_finer(index_y))
if (slope<0) then 
    q_menu_fun = value_left
else
    q_menu_fun = max(min(value_left + slope * (y - y_grid_finer(index_y)), EXP(-r)),0d+0)
end if
!q_menu_fun = interpolate_finer(b_local,y,q_menu_matrix)
end        
    
DOUBLE PRECISION function EV_excl_fun(b, y)
USE param
DOUBLE PRECISION, INTENT(IN) :: b, y
DOUBLE PRECISION ::  slope, DCSVAL, value_left, value_right, b_local, interpolate_finer
INTEGER :: index_y, NINTV
EXTERNAL DCSVAL, interpolate_finer

!!RESTRICT TO b INSIDE THE GRID SPACE. OTHERWISE, INTERPOLATION PROCEDURE FAILS
b_local = min(max(b, b_inf+1d-12), b_sup-1d-12)
NINTV = b_num_finer - 1
! THIS IS FOR SPLINE
index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))
value_left  = DCSVAL(b, NINTV, break_matrix_ev_excl(:, index_y), coeff_matrix_ev_excl(:,:, index_y))
value_right = DCSVAL(b, NINTV, break_matrix_ev_excl(:, index_y+1), coeff_matrix_ev_excl(:,:, index_y+1))
slope = (value_right - value_left)/ (y_grid_finer(index_y+1) - y_grid_finer(index_y))
EV_excl_fun = value_left + slope * (y - y_grid_finer(index_y))
!EV_excl_fun = interpolate_finer(b_local,y,EV_excl_matrix)
end  
    

DOUBLE PRECISION function EW_excl_fun(b, y)
USE param
DOUBLE PRECISION, INTENT(IN) :: b, y
DOUBLE PRECISION ::  slope, DCSVAL, value_left, value_right, b_local, interpolate_finer
INTEGER :: index_y, NINTV
EXTERNAL DCSVAL, interpolate_finer

!!RESTRICT TO b INSIDE THE GRID SPACE. OTHERWISE, INTERPOLATION PROCEDURE FAILS
b_local = min(max(b, b_inf+1d-12), b_sup-1d-12)
NINTV = b_num_finer - 1
! THIS IS FOR SPLINE
index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))
value_left  = DCSVAL(b, NINTV, break_matrix_ew_excl(:, index_y), coeff_matrix_ew_excl(:,:, index_y))
value_right = DCSVAL(b, NINTV, break_matrix_ew_excl(:, index_y+1), coeff_matrix_ew_excl(:,:, index_y+1))
slope = (value_right - value_left)/ (y_grid_finer(index_y+1) - y_grid_finer(index_y))
EW_excl_fun = value_left + slope * (y - y_grid_finer(index_y))
!EW_excl_fun = interpolate_finer(b_local,y,EW_excl_matrix)
end      
    
    
DOUBLE PRECISION function EV_fun(b, y)
USE param
DOUBLE PRECISION, INTENT(IN) :: b, y
DOUBLE PRECISION ::  slope, DCSVAL, value_left, value_right, b_local, interpolate_finer
INTEGER :: index_y, NINTV
EXTERNAL DCSVAL, interpolate_finer

!!RESTRICT TO b INSIDE THE GRID SPACE. OTHERWISE, INTERPOLATION PROCEDURE FAILS
b_local = min(max(b, b_inf+1d-12), b_sup-1d-12)
NINTV = b_num_finer - 1
! THIS IS FOR SPLINE
index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))
value_left  = DCSVAL(b, NINTV, break_matrix_ev(:, index_y), coeff_matrix_ev(:,:, index_y))
value_right = DCSVAL(b, NINTV, break_matrix_ev(:, index_y+1), coeff_matrix_ev(:,:, index_y+1))
slope = (value_right - value_left)/ (y_grid_finer(index_y+1) - y_grid_finer(index_y))
EV_fun = value_left + slope * (y - y_grid_finer(index_y))
!EV_fun = interpolate_finer(b_local,y,EV_matrix)
end
    
DOUBLE PRECISION function EW_fun(b, y)
USE param
DOUBLE PRECISION, INTENT(IN) :: b, y
DOUBLE PRECISION ::  slope, DCSVAL, value_left, value_right, b_local, interpolate_finer
INTEGER :: index_y, NINTV
EXTERNAL DCSVAL, interpolate_finer

!!RESTRICT TO b INSIDE THE GRID SPACE. OTHERWISE, INTERPOLATION PROCEDURE FAILS
b_local = min(max(b, b_inf+1d-12), b_sup-1d-12)
NINTV = b_num_finer - 1
! THIS IS FOR SPLINE
index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))
value_left  = DCSVAL(b, NINTV, break_matrix_ew(:, index_y), coeff_matrix_ew(:,:, index_y))
value_right = DCSVAL(b, NINTV, break_matrix_ew(:, index_y+1), coeff_matrix_ew(:,:, index_y+1))
slope = (value_right - value_left)/ (y_grid_finer(index_y+1) - y_grid_finer(index_y))
EW_fun = value_left + slope * (y - y_grid_finer(index_y))
!EW_fun = interpolate_finer(b_local,y,EW_matrix)
end     
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIN NUEVO!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!COMPUTE MATRICES AND SPLINE COEFFICIENTS FOR q(b, y, s) and E[V(b, y') | y)]
subroutine compute_q_ev
USE param
INTEGER :: i_b, i_y, i, ILEFT, IRIGHT, LDF, KXORD, KYORD
parameter (KXORD =3, KYORD = 3)
DOUBLE PRECISION :: b_long, y_threshold, cdf_threshold, y_fun, DNORDF, exp_v_next, scalar, y_min, y_max, DNORIN, exp_v_excl_ends, exp_v_excl,&
                    acum1, acum2, acum3, cdf, cdf1, y_next, y_next1, value_next, v0_fun, v1_fun, q_fun, b_long1, b_next_post_def, acum,&
                    interpolate3, exp_v_excl1, exp_v_excl_ends1, w0_fun, w1_fun, acum1_w, acum2_w, value_next_w,&
                    DCSVAL, FDATA(b_num_finer), DLEFT, DRIGHT, BREAK_GRID(b_num_finer), CSCOEF(4,b_num_finer)
EXTERNAL y_fun, DNORDF, v0_fun, v1_fun, DCSDEC, DNORIN, q_fun, DBS2IN, interpolate3, DCSVAL, w0_fun, w1_fun

do i_y = 1,y_num_finer
    !WRITE(nout, *) i_y
    y_initial = y_grid_finer(i_y)
    do i_b = 1, b_num_finer
        b_long = b_grid_finer(i_b)
        !COMPUTE EXPECTED VALUE FUNCTION FOR THE NEXT PERIOD IF IT CHOOSES A PORTFOLIO (b_long)
        !WHEN CURRENT INCOME = y_initial
        y_threshold   = y_fun(b_long)
        cdf_threshold = DNORDF((y_threshold - rho*y_initial - (1-rho)* mean_y)/std_eps)
        exp_v_next = 0d+0
        acum1=0d+0
        acum2=0d+0
        acum1_w=0d+0
        acum2_w=0d+0
        if (cdf_threshold < cdf_inf) then
            do i=1,quad_num_v
                cdf = 0.5*(quad_x_v(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
                y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                value_next = prob_office*v0_fun(b_long, y_next) + (1-prob_office)*w0_fun(b_long, y_next)
                acum1 = acum1 + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next
                value_next_w = prob_office*w0_fun(b_long, y_next) + (1-prob_office)*v0_fun(b_long, y_next)
                acum1_w = acum1_w + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next_w
            end do
        elseif(cdf_threshold > cdf_sup) then
            do i=1,eps_num
                cdf = 0.5*(quad_x_v(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
                y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                value_next = prob_office*v1_fun(b_long, y_next) + (1-prob_office)*w1_fun(b_long, y_next)
                acum1 = acum1 + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next
                value_next_w = prob_office*w1_fun(b_long, y_next) + (1-prob_office)*v1_fun(b_long, y_next)
                acum1_w = acum1_w + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next_w
            end do
        else
            do i=1,eps_num
                cdf = 0.5*(quad_x_v(i) + 1)*(cdf_threshold - cdf_inf) + cdf_inf
                y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                value_next = prob_office*v1_fun(b_long, y_next) + (1-prob_office)*w1_fun(b_long, y_next)
                acum1 = acum1 + 0.5*(cdf_threshold - cdf_inf)*quad_w_v(i) * value_next
                value_next_w = prob_office*w1_fun(b_long, y_next) + (1-prob_office)*v1_fun(b_long, y_next)
                acum1_w = acum1_w + 0.5*(cdf_threshold - cdf_inf)*quad_w_v(i) * value_next_w
        
                !compute integral for y > y_threshold
                cdf1 = 0.5*(quad_x_v(i) + 1)*(cdf_sup - cdf_threshold) + cdf_threshold
                y_next1 = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf1) * std_eps !DCSVAL(cdf1, NINTV, BREAK_eps, CSCOEF_eps)
                value_next = prob_office*v0_fun(b_long, y_next1) + (1-prob_office)*w0_fun(b_long, y_next1)
                acum2 = acum2 + 0.5*(cdf_sup - cdf_threshold)*quad_w_v(i) * value_next
                value_next_w = prob_office*w0_fun(b_long, y_next1) + (1-prob_office)*v0_fun(b_long, y_next1)
                acum2_w = acum2_w + 0.5*(cdf_sup - cdf_threshold)*quad_w_v(i) * value_next_w
                !WRITE(119, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') y_next, v1_fun(y_next), y_next1, value_next
            end do
        end if
        acum = (acum1+acum2)/(cdf_sup - cdf_inf)
        EV_matrix(i_b, i_y) = acum
        EW_matrix(i_b, i_y) = (acum1_w+acum2_w)/(cdf_sup - cdf_inf)
        q_menu_matrix(i_b, i_y) = q_fun(b_long, y_initial)
        !write(nout, '(G15.6, X, G15.6, X, G15.6, X, G15.6, X, G15.6)') b_long, y_initial, acum, v0_fun(b_long, y_initial), v1_fun(b_long, y_initial)
        !if (i_b == 40) then
        !    pause
        !    write(nout, '(G15.6, X, G15.6, X, G15.6, X, G15.6, X, G15.6)') b_long, y_initial, acum, v0_fun(b_long, y_initial), v1_fun(b_long, y_initial)
        !end if
    end do
    acum =0d+0
    acum1=0d+0
    acum2=0d+0
    acum1_w=0d+0
    acum2_w=0d+0
    b_long  = 0d+0
    do i=1,quad_num_v
        !CONTINUES EXCLUDED IN THE NEXT PERIOD (NO BORROWING AND DEFAULT COST)
        cdf = 0.5*(quad_x_v(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
        y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
        value_next = prob_office*v0_fun(b_long, y_next) + (1-prob_office)*w0_fun(b_long, y_next)
        acum1 = acum1 + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next
        value_next_w = prob_office*w0_fun(b_long, y_next) + (1-prob_office)*v0_fun(b_long, y_next)
        acum1_w = acum1_w + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next_w
        
        !EXCLUSION ENDS IN THE NEXT PERIOD (STARTS WITH DEBT AND NO DEFAULT COST)
        value_next = prob_office*v1_fun(b_long, y_next) + (1-prob_office)*w1_fun(b_long, y_next)
        acum2 = acum2 + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next
        value_next_w = prob_office*w1_fun(b_long, y_next) + (1-prob_office)*v1_fun(b_long, y_next)
        acum2_w = acum2_w + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next_w
    end do
    acum = ((1d+0 - prob_excl_end)*acum2 + prob_excl_end*acum1)/(cdf_sup - cdf_inf)
    EV_excl_matrix(:, i_y) = acum
    EW_excl_matrix(:, i_y) = ((1d+0 - prob_excl_end)*acum2_w + prob_excl_end*acum1_w)/(cdf_sup - cdf_inf)
end do        

do i_y = 1,y_num_finer
    FDATA = EV_excl_matrix(:,i_y)
    ILEFT  = 0
    IRIGHT = 0
    DLEFT = (FDATA(2)-FDATA(1)) / (b_grid_finer(2) - b_grid_finer(1))
    DRIGHT = (FDATA(b_num_finer)-FDATA(b_num_finer-1)) / (b_grid_finer(b_num_finer) - b_grid_finer(b_num_finer-1))
    call  DCSDEC (b_num_finer, b_grid_finer, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
    break_matrix_EV_excl(:, i_y) = BREAK_GRID
    coeff_matrix_EV_excl(:,:, i_y) = CSCOEF
        
    FDATA = EV_matrix(:,i_y)
    ILEFT  = 0
    IRIGHT = 0
    DLEFT = (FDATA(2)-FDATA(1)) / (b_grid_finer(2) - b_grid_finer(1))
    DRIGHT = (FDATA(b_num_finer)-FDATA(b_num_finer-1)) / (b_grid_finer(b_num_finer) - b_grid_finer(b_num_finer-1))
    call  DCSDEC (b_num_finer, b_grid_finer, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
    break_matrix_EV(:, i_y) = BREAK_GRID
    coeff_matrix_EV(:,:, i_y) = CSCOEF
    
    FDATA = EW_excl_matrix(:,i_y)
    ILEFT  = 0
    IRIGHT = 0
    DLEFT = (FDATA(2)-FDATA(1)) / (b_grid_finer(2) - b_grid_finer(1))
    DRIGHT = (FDATA(b_num_finer)-FDATA(b_num_finer-1)) / (b_grid_finer(b_num_finer) - b_grid_finer(b_num_finer-1))
    call  DCSDEC (b_num_finer, b_grid_finer, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
    break_matrix_EW_excl(:, i_y) = BREAK_GRID
    coeff_matrix_EW_excl(:,:, i_y) = CSCOEF
        
    FDATA = EW_matrix(:,i_y)
    ILEFT  = 0
    IRIGHT = 0
    DLEFT = (FDATA(2)-FDATA(1)) / (b_grid_finer(2) - b_grid_finer(1))
    DRIGHT = (FDATA(b_num_finer)-FDATA(b_num_finer-1)) / (b_grid_finer(b_num_finer) - b_grid_finer(b_num_finer-1))
    call  DCSDEC (b_num_finer, b_grid_finer, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
    break_matrix_EW(:, i_y) = BREAK_GRID
    coeff_matrix_EW(:,:, i_y) = CSCOEF
    
    FDATA = q_menu_matrix(:,i_y)
    ILEFT  = 0
    IRIGHT = 0
    DLEFT = (FDATA(2)-FDATA(1)) / (b_grid_finer(2) - b_grid_finer(1))
    DRIGHT = (FDATA(b_num_finer)-FDATA(b_num_finer-1)) / (b_grid_finer(b_num_finer) - b_grid_finer(b_num_finer-1))
    call  DCSDEC (b_num_finer, b_grid_finer, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
    break_matrix_q_menu(:, i_y) = BREAK_GRID
    coeff_matrix_q_menu(:,:, i_y) = CSCOEF
end do

!if (DMOD(counter, 5d+0)==zero) then
!open(100, FILE ='graphs\coeff_excl.txt',STATUS='replace')
!open(101, FILE ='graphs\coeff_EV.txt',STATUS='replace')
!open(102, FILE ='graphs\coeff_q.txt',STATUS='replace')
!
!do i_access = 1, 2
!    do i_y = 1,y_num_finer
!        do i_b_short = 1,b_num_short_finer
!            WRITE(100, '(F15.6, X, F15.6, X, F15.6, X, F15.6, X, F15.6)') break_matrix_EV_excl(i_access, i_y, i_b_short), &
!            coeff_matrix_EV_excl(i_access, i_y, 1, i_b_short), coeff_matrix_EV_excl(i_access, i_y, 2, i_b_short), &
!            coeff_matrix_EV_excl(i_access, i_y, 3, i_b_short), coeff_matrix_EV_excl(i_access, i_y, 4, i_b_short)
!        end do
!        
!        do i_b_long = 1,b_num_long_finer
!            do i_b_short = 1,b_num_short_finer
!                WRITE(101, '(F15.6, X, F15.6, X, F15.6, X, F15.6, X, F15.6)') break_matrix_EV(i_access, i_y, i_b_long, i_b_short), &
!                coeff_matrix_EV(i_access, i_y, i_b_long, 1, i_b_short), coeff_matrix_EV(i_access, i_y, i_b_long, 2, i_b_short), &
!                coeff_matrix_EV(i_access, i_y, i_b_long, 3, i_b_short), coeff_matrix_EV(i_access, i_y, i_b_long, 4, i_b_short)
!                
!                WRITE(102, '(F15.6, X, F15.6, X, F15.6, X, F15.6, X, F15.6)') break_matrix_q_menu(i_access, i_y, i_b_long, i_b_short), &
!                coeff_matrix_q_menu(i_access, i_y, i_b_long, 1, i_b_short), coeff_matrix_q_menu(i_access, i_y, i_b_long, 2, i_b_short), &
!                coeff_matrix_q_menu(i_access, i_y, i_b_long, 3, i_b_short), coeff_matrix_q_menu(i_access, i_y, i_b_long, 4, i_b_short)
!            end do
!        end do
!    end do
!end do
!CLOSE(100)
!CLOSE(101)
!CLOSE(102)
!end if
end subroutine    


subroutine iterate
USE param
INTEGER :: d
DOUBLE PRECISION :: y_valor, b_valor, b_valor_def, v_valor, convergence, criteria, deviation, q_fun, b_next, g,&
                    b0_next, b1_next, b, w_value, v_valor_excl, objective_function, dev_q_paid, dev_q_scheme, dev_v,&
                    v0_value, v1_value, q, interpolate, b_next1, b_next0,  b_grid_local(b_num), b_inf_local, b_inf_old,&
                    FDATA(b_num), DLEFT, DRIGHT, BREAK_GRID(b_num), CSCOEF(4,b_num), b_tirar, DCSVAL, &
                    epsilon_grid(cdf_num), dev_q_nodef, objective_excl, b_adj,  dev_q_paid_free, objective_excl_w, objective_function_w

DOUBLE PRECISION, DIMENSION(b_num, y_num) :: v0_matrix_new, v1_matrix_new, q_matrix_new, default_decision_new,&
                                             q_scheme, q_scheme_new, v_excl_matrix_new, w0_matrix_new, w1_matrix_new, w_excl_matrix_new
DOUBLE PRECISION, DIMENSION(b_num, y_num) :: v_matrix_new, dev_matrix, q_matrix_nodef_new, w_matrix_new
INTEGER i_b, i_y, i, j, i_def_opt, i_b_zero, indices(1:2), indice_b, num_tirar, ILEFT, IRIGHT,NINTV,&
        b_graph_num, i_m
EXTERNAL q_fun, objective_function, DCSDEC, DCSVAL, DNORIN, objective_excl, b_adj, objective_excl_w, objective_function_w


!compute spline coefficients of inverse (gaussian cdf)
!COMPUTE SPLINE COEFFICIENT OF INVERSE cdf(epsilon)
!USED LATER TO COMPUTE EXPECTATIONS.
do i=1,cdf_num
   epsilon_grid(i) = DNORIN(cdf_grid(i)) * std_eps
end do

NINTV  = cdf_num - 1
ILEFT  = 0
IRIGHT = 0
DLEFT = (epsilon_grid(2) - epsilon_grid(1)) / (cdf_grid(2) - cdf_grid(1))
DRIGHT = (epsilon_grid(cdf_num) - epsilon_grid(cdf_num-1)) / (cdf_grid(cdf_num) - cdf_grid(cdf_num-1))
call  DCSDEC (cdf_num, cdf_grid, epsilon_grid, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_eps, CSCOEF_eps)

criteria = 1d-4   !CRITERIA FOR CONVERGENCE
convergence = -1
ERRREL = 1d-10    !PRECISION WITH WHICH THE SPLINE COEFFICIENTS ARE COMPUTED

b_inf_old = b_inf
v1_matrix_new = v1_matrix
v0_matrix_new = v0_matrix
v_matrix_new = v_matrix

do WHILE(convergence<0)
    deviation = 0d+0
    counter = counter + 1
    dev_v = 0d+0
    dev_q_scheme = 0d+0
    dev_q_paid   = 0d+0
    dev_q_nodef  = 0d+0
    dev_q_nodil  = 0d+0
    dev_q_paid_free = 0d+0
    do i_y = 1,y_num
        FDATA = v0_matrix(:,i_y)
	    ILEFT  = 0
	    IRIGHT = 0
	    DLEFT = (FDATA(2)-FDATA(1)) / (b_grid(2) - b_grid(1))
	    DRIGHT = (FDATA(b_num)-FDATA(b_num-1)) / (b_grid(b_num) - b_grid(b_num-1))
	    call  DCSDEC (b_num, b_grid, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
        break_matrix(i_y, :) = BREAK_GRID
        coeff_matrix(i_y, :,:) = CSCOEF
        
        FDATA = v1_matrix(:,i_y)
	    ILEFT  = 0
	    IRIGHT = 0
	    DLEFT = (FDATA(2)-FDATA(1)) / (b_grid(2) - b_grid(1))
	    DRIGHT = (FDATA(b_num)-FDATA(b_num-1)) / (b_grid(b_num) - b_grid(b_num-1))
	    call  DCSDEC (b_num, b_grid, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
        break_matrix_v1(i_y, :) = BREAK_GRID
        coeff_matrix_v1(i_y, :,:) = CSCOEF
        
     !   FDATA = w_matrix(:,i_y)
	    !ILEFT  = 0
	    !IRIGHT = 0
	    !DLEFT = (FDATA(2)-FDATA(1)) / (b_grid(2) - b_grid(1))
	    !DRIGHT = (FDATA(b_num)-FDATA(b_num-1)) / (b_grid(b_num) - b_grid(b_num-1))
	    !call  DCSDEC (b_num, b_grid, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
     !   break_matrix_w(i_y, :) = BREAK_GRID
     !   coeff_matrix_w(i_y, :,:) = CSCOEF
        
        FDATA = w0_matrix(:,i_y)
	    ILEFT  = 0
	    IRIGHT = 0
	    DLEFT = (FDATA(2)-FDATA(1)) / (b_grid(2) - b_grid(1))
	    DRIGHT = (FDATA(b_num)-FDATA(b_num-1)) / (b_grid(b_num) - b_grid(b_num-1))
	    call  DCSDEC (b_num, b_grid, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
        break_matrix_w0(i_y, :) = BREAK_GRID
        coeff_matrix_w0(i_y, :,:) = CSCOEF
        
        FDATA = w1_matrix(:,i_y)
	    ILEFT  = 0
	    IRIGHT = 0
	    DLEFT = (FDATA(2)-FDATA(1)) / (b_grid(2) - b_grid(1))
	    DRIGHT = (FDATA(b_num)-FDATA(b_num-1)) / (b_grid(b_num) - b_grid(b_num-1))
	    call  DCSDEC (b_num, b_grid, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
        break_matrix_w1(i_y, :) = BREAK_GRID
        coeff_matrix_w1(i_y, :,:) = CSCOEF
        
        FDATA = q_matrix_nodef(:,i_y)
        ILEFT  = 0
        IRIGHT = 0
        call  DCSDEC (b_num, b_grid, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
        break_matrix_q(i_y, :) = BREAK_GRID
        coeff_matrix_q(i_y, :,:) = CSCOEF
    end do
    call compute_q_ev


    WRITE(nout, *) 'Finished'
    do i_y = 1,y_num
        !WRITE(nout, *) i_y
        i_y_global = i_y
        y_initial = y_grid(i_y)
        b_initial = 0d+0
        i_default_global=2 !Country defaults
        v_valor_excl = - objective_excl(b_initial)
        !indicator_tirar = 1
        !call optimize(b_next1, v_valor)
        !indicator_tirar = 0
        !WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8)') b_next1, q_fun(b_next1, y_initial), v_valor
        !pause
        v1_matrix_new(:, i_y) = v_valor_excl
        v_excl_matrix_new(:, i_y) = v_valor_excl
        w1_matrix_new(:, i_y) = - objective_excl_w(b_initial)
        !WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8)') objective_excl_w(b_initial), v_valor_excl
        w_excl_matrix_new(:, i_y) = - objective_excl_w(b_initial)
        do i_b = 1, b_num !1,i_b_zero
            b_initial = b_grid(i_b)
            i_default_global=1 !Country does not default
            call optimize(b_next0, v_valor)
            v0_matrix_new(i_b, i_y) = v_valor
            w0_matrix_new(i_b, i_y) = - objective_function_w(b_next0)
            q_matrix_nodef_new(i_b, i_y) = q_fun(b_next0, y_initial)
            indicator_tirar = 0
            if (v1_matrix_new(i_b, i_y) > v0_matrix_new(i_b, i_y)) THEN !IT IS OPTIMAL TO DEFAULT
                b_next_matrix(i_b, i_y) = zero
                default_decision_new(i_b, i_y) = 2
                v_matrix_new(i_b, i_y) = v1_matrix_new(i_b, i_y)
                q_matrix_new(i_b, i_y) = 0
                w_matrix_new(i_b, i_y) = w1_matrix_new(i_b, i_y)
            else
                b_next_matrix(i_b, i_y) = b_next0
                default_decision_new(i_b, i_y) = 1
                v_matrix_new(i_b, i_y) = v0_matrix_new(i_b, i_y)
                w_matrix_new(i_b, i_y) = w0_matrix_new(i_b, i_y)
                q_matrix_new(i_b, i_y) = q_fun(b_next_matrix(i_b, i_y), y_initial)
            END if
            !UPDATE MATRIX OF PAID q AT EACH STATE
            !WRITE(nout, '(I3, X, I3, X, F10.6, X, F10.6, X, F10.6, X, F10.6)') i_y, i_b,&
            !q_matrix(i_b, i_y) , q_matrix_new(i_b, i_y),q_matrix_new(i_b, i_y) - q_matrix(i_b, i_y), b_next_matrix(i_b,i_y)
            dev_v =        max(ABS(v_matrix_new(i_b, i_y) - v_matrix(i_b, i_y)), dev_v)
            dev_q_paid =   MAX(ABS(q_matrix_new(i_b, i_y) - q_matrix(i_b, i_y)), dev_q_paid)
            !deviation = MAX(deviation, MAX(dev_v, dev_q_paid, dev_q_paid_free))
            deviation = MAX(deviation, dev_v)
            dev_matrix(i_b, i_y) = (q_matrix_new( i_b, i_y) - q_matrix( i_b, i_y))
            !write(nout, *) i_b, i_y, dev_v, dev_q_paid
          end do
    end do
!indices = MAXLOC(dev_matrix)
! CLOSE(117)
    WRITE(nout, '(F12.8, X, F4.0, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F10.6)') deviation, counter, dev_v, dev_q_paid

!3) SAVE RESULTS OF THE CURRENT ITERATION

    open (10, FILE='graphs\v.txt',STATUS='replace')
    open (11, FILE='graphs\default.txt',STATUS='replace')
    !open (12, FILE='graphs\q.txt',POSITION='append')
    open (12, FILE='graphs\q.txt', STATUS='replace')
    open (13, FILE='graphs\b_next.txt',STATUS='replace')
    open (14, FILE='graphs\dev.txt',STATUS='replace')
    !open (16, FILE='graphs\q_paid.txt',POSITION='append')
    open (16, FILE='graphs\q_paid.txt', STATUS='replace')
    open (17, FILE='graphs\counter.txt', STATUS='replace')
    open (18, FILE='graphs\w.txt',STATUS='replace')
    open (110, FILE='graphs\b_grid.txt',STATUS='replace')
    open (113, FILE='graphs\bounds.txt',STATUS='replace')
    WRITE(113, '(F15.11, X, F15.11)') b_inf, b_sup
    WRITE(113, '(F15.11, X, F15.11)') y_inf, y_sup
    WRITE(17, '(F10.1)') counter
    do i_b = 1,b_num
        WRITE(110, '(F12.8)') b_grid(i_b)
        do i_y = 1,y_num
            i_y_global = i_y
            y_initial = y_grid(i_y)
            b_initial = b_grid(i_b)
            i_b_global = i_b

            indicator_tirar=0
            WRITE(10, '(F15.10, X, F15.10, X, F15.10, X, F15.10)') v_matrix_new(i_b, i_y), &
            v0_matrix_new(i_b, i_y), v1_matrix_new(i_b, i_y), v_excl_matrix_new(i_b, i_y)
            WRITE(11, '(F6.2)') default_grid(default_decision_new(i_b, i_y))
            d = default_decision_new(i_b, i_y)
            b = b_grid(i_b)
            g = y_grid(i_y)
            WRITE(12, '(F15.11)') q_fun(b_grid(i_b), y_grid(i_y))
            WRITE(13, '(F15.11)') b_next_matrix(i_b, i_y)
            WRITE(14, '(F15.11)') dev_matrix(i_b, i_y)
            WRITE(16, '(F15.11, X, F15.11)') q_matrix_new(i_b, i_y), q_matrix_nodef_new(i_b, i_y)
            WRITE(18, '(F15.10, X, F15.10, X, F15.10, X, F15.10)') w_matrix_new(i_b, i_y), &
            w0_matrix_new(i_b, i_y), w1_matrix_new(i_b, i_y), w_excl_matrix_new(i_b, i_y)
            indicator_tirar=0
        end do
    end do
    CLOSE(10)
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
    CLOSE(14)
    CLOSE(16)
    CLOSE(110)
    CLOSE(113)
    CLOSE(17)
    CLOSE(18)
    !pause
    !UPDATE VALUES OF MATRICES

    default_decision = default_decision_new
    q_matrix = q_matrix_new
    q_matrix_nodef = q_matrix_nodef_new
    v0_matrix = v0_matrix_new
    v1_matrix = v1_matrix_new
    v_excl_matrix = v_excl_matrix_new
    v_matrix = v_matrix_new
    w0_matrix = w0_matrix_new
    w1_matrix = w1_matrix_new
    w_excl_matrix = w_excl_matrix_new
    w_matrix = w_matrix_new
    !PRINT*,'deviation =', deviation

    if (deviation < criteria) then
        convergence =1     
        !call b_adj
        call simulate
    end if
end do
!FINALLY, STORE VALUE FUNCTIONS, POLICY FUNCTIONS AND PRICES USING A FINER GRID.
!THIS HELPS TO VISUALIZE THE RESULTS.
end subroutine
    

   
subroutine simulate
USE param
integer :: period_num, i,j,k, gov_type, random_num, ivalue, sample_num, MAXFN
parameter (sample_num = 250, period_num=801)
DOUBLE PRECISION :: random_matrix(period_num, sample_num, 3), random_vector(1:3*sample_num*period_num), &
                    z(period_num,sample_num), b(period_num+1,sample_num), b_adj_fun, b_adjust(period_num+1,sample_num),&
                    q(period_num,sample_num), eps, b_interp(b_num, y_num), v_def, v_no_def, b_next, q_fun, &
                    q_interp(b_num, y_num),q_paid, b_zero, c(period_num,sample_num), tb(period_num,sample_num),&
                    y(period_num,sample_num), STEP, BOUND, XACC, objective_function, b_next_guess, gamma,&
                    v_valor, v0_fun, v1_fun, DNORIN
INTEGER :: excl(period_num, sample_num), d(period_num, sample_num)
EXTERNAL RNSET, DRNUN, DNORIN, q_fun, DUVMIF, objective_function, v0_fun, v1_fun

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
!open (UNIT=1, FILE="graphs\r_vector.txt", status = 'replace')
!   do i=1,period_num*sample_num*2
!      write(1,'(F12.8)') random_vector(i)
!      !READ
!      end do
!CLOSE(1)

!First column of random_matrix is used to generate transitory shocks.
!Second column of random_matrix is used to generate type changes.
do j=1,sample_num
    do i=1,period_num
       random_matrix(i,j,1)=random_vector((j-1)*period_num+i)
       random_matrix(i,j,2)=random_vector(period_num*sample_num + (j-1)*period_num+i)
       random_matrix(i,j,3)=random_vector(2*period_num*sample_num + (j-1)*period_num+i)
    end do
end do

i_type_global =1
    !Retrieve value functions and policy functions for a type t gov't that has made a default d(i-1) decision
    !in the previous period

open (UNIT=21, FILE="graphs\data_sim.txt", status = 'replace')
open (UNIT=22, FILE="graphs\def_per.txt", status = 'replace')
open (UNIT=23, FILE="graphs\param.txt", status = 'replace')
open (UNIT=24, FILE="graphs\tirar.txt", status = 'replace')
open (UNIT=25, FILE="graphs\delta.txt", status = 'replace')

WRITE(23, '(I10)') period_num-1
WRITE(23, '(I10)') sample_num
CLOSE(23)

WRITE(25, '(F15.10)') delta
CLOSE(25)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
i_default_global = 1  !USED TO INVOKE WHICH CHEBYCHEV MATRIX TO USE
                      !NEED TO MODIFY THIS WHEN DEFAULT AFFECTS OUTPUT REGARDLESS OF THE EXCLUSION STATUS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=1,sample_num   !SOLVE FOR SAMPLE j
    WRITE(nout, *) j
    !Set initial values
    z(1,j) = mean_y
    y(1,j) = EXP(z(1,j))
    b(1,j) = zero
    b(2,j) = zero
    d(1,j) = 1   !NO DEFAULT IN FIRST PERIOD
    excl(1,j) = 0 !Country is not excluded in the first period
    do i=2,period_num
        !epsilon = realization of standard gaussian * standard deviation
        eps = DNORIN(random_matrix(i,j,1)) * std_eps !GROWTH SHOCK
        z(i,j) = rho*z(i-1,j) + (1-rho)*mean_y + eps   !current output = ex ante mean + epsilon
        b_initial = b(i,j)
        y_initial = z(i,j)
        v_no_def = v0_fun(b_initial, y_initial)
        v_def = v1_fun(b_initial, y_initial)
        !WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8)') b_initial, y_initial, v_no_def, v_def
        if (excl(i-1,j) < 0.5) then !Country was not excluded in the previous period
            if (v_def>v_no_def) then  !COUNTRY DEFAULTS
                d(i,j) = 2
                i_default_global = 2
                excl(i,j) = 1
                b(i+1,j) = 0d+0
                !b_global = b_next
                q(i,j) = 1d+0 / (delta + r) !q_fun(b_next, z(i,j))
                y(i,j) = EXP(z(i,j)) * (1-d0) - d1*EXP(z(i,j))**2d+0
                c(i,j) = y(i,j)! - q(i,j)* b(i+1,j) 
                tb(i,j) = y(i,j) - c(i,j)
                WRITE(22, '(I7)') i-1 !NEED TO SUBSTRACT 1. REASON: files start saving data on period 2
            else
                d(i,j) = 1   !COUNTRY DOES NOT DEFAULT
                excl(i,j) = 0
                i_default_global = 1
                call optimize(b_next, v_valor)
                b(i+1,j) = b_next
                !Compute the bond price paid when an amount b_next is issued.
                !the bond price depends on the current default decision and output (both will affect output tomorrow)
                !It also depends on the current type in power.
                q(i,j) = q_fun(b_next, z(i,j))
                y(i,j) = EXP(z(i,j))
                c(i,j) = y(i,j) + coupon*b(i,j) - q(i,j)*(b(i+1,j) - b(i,j)*(1d+0 - delta))
                tb(i,j) = y(i,j) - c(i,j)
            end if
        else !Country was excluded in previous period
            if (random_matrix(i,j,2) <= prob_excl_end) then !exclusion ends today
                d(i,j) = 1 !country does not default
                excl(i,j) = 0
                i_default_global = 1
                call optimize(b_next, v_valor)
                b(i+1,j) = b_next
                q(i,j) = q_fun(b_next, z(i,j))
                y(i,j) = EXP(z(i,j))
                c(i,j) = y(i,j) + coupon*b(i,j) - q(i,j)*(b(i+1,j) - b(i,j)*(1d+0 - delta))
                tb(i,j) = y(i,j) - c(i,j)
            else !exclusion continues
                d(i,j) = 1
                excl(i,j) = 1
                i_default_global = 1
                b(i+1,j) = 0d+0
                !b_global = b_next
                q(i,j) = 1d+0 / (delta + r) !q_fun(b_next, z(i,j))
                y(i,j) = EXP(z(i,j)) * (1-d0) - d1*EXP(z(i,j))**2d+0  !output loss while exclusion
                c(i,j) = y(i,j)! - q(i,j)* b(i+1,j)
                tb(i,j) = y(i,j) - c(i,j)
            end if
        end if
        WRITE(21, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, I3, X, I3)') &
        LOG(y(i,j)), b(i,j), q(i,j), LOG(c(i,j)), tb(i,j)/y(i,j), d(i,j), excl(i,j)
        !WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, I3, X, I3, X, F8.5)') &
        !LOG(y(i,j)), b(i+1,j), q(i,j), LOG(c(i,j)), tb(i,j)/y(i,j), d(i,j)
    end do
end do
CLOSE(21)
CLOSE(22)
CLOSE(24)
end subroutine
    
program main
include 'link_fnl_static.h'
USE param
DOUBLE PRECISION :: y_valor, b_valor, f_valor, q_fun, u_fun, start_time, end_time, indicator_external, def, u_fun_w
INTEGER  i_b, i_y
EXTERNAL q_fun, u_fun, u_fun_w


call cpu_time(start_time)
call quadrature

indicator_external = 1 !FROM EXTERNAL FILE       ! THIS IS TO TAKE THE INITIAL GUESS FROM AN OUTSIDE FILE ( FROM PREVIOUS ITERATIONS )
if (indicator_external < 0.5) then
    !INITIAL BOUNDS ON GRID FOR b
    call compute_grid
    counter = 0
    open (10, FILE='graphs\v.txt',STATUS='replace')
    open (11, FILE='graphs\default.txt',STATUS='replace')
    do i_b = 1,b_num
        b_initial = b_grid(i_b)
        !WRITE(nout, '(A10, X, A10, X, A10, X, A10)') 'y', 'b', 'c1', 'c0'
        do i_y = 1,y_num
            y_initial = y_grid(i_y)
            v1_matrix(i_b, i_y) = u_fun(MIN(lambda, exp(y_initial)))
            w1_matrix(i_b, i_y) = u_fun_w(MIN(lambda, exp(y_initial)))
            !v1_matrix(i_b, i_y) = u_fun(EXP(y_grid(i_y)) + b_grid(i_b))
            v0_matrix(i_b, i_y) = u_fun(EXP(y_grid(i_y)) + coupon* b_grid(i_b)*(1-delta))
            v_matrix(i_b, i_y) = MAX(v1_matrix(i_b, i_y), v0_matrix(i_b,i_y))
            w0_matrix(i_b, i_y) = u_fun_w(EXP(y_grid(i_y)) + coupon* b_grid(i_b)*(1-delta))
            w_matrix(i_b, i_y) = MAX(w1_matrix(i_b, i_y), w0_matrix(i_b,i_y))
            if (v1_matrix(i_b, i_y) <= v0_matrix(i_b, i_y)) then
            default_decision(i_b, i_y) = 1
            else
            default_decision(i_b, i_y) = 2
            end if
            WRITE(10, '(F15.10, X, F15.10, X, F15.10, X, F15.10)') v_matrix(i_b, i_y), &
                v0_matrix(i_b, i_y), v1_matrix(i_b, i_y), v_excl_matrix(i_b, i_y)            
            WRITE(11, '(F6.2)') default_grid(default_decision(i_b, i_y))
            q_matrix(i_b, i_y) = 0
            q_matrix_nodef(i_b, i_y) = 0
        end do
    end do
    CLOSE(10)
    CLOSE(11)
else !READ DATA FROM EXTERNAL FILES
    open (10, FILE='graphs\v.txt')
    !open (12, FILE='graphs\default.txt')
    open (16, FILE='graphs\q_paid.txt')
    open (13, FILE='graphs\b_next.txt')
    open (17, FILE='graphs\counter.txt')
    READ(17, '(F10.1)') counter
    call compute_grid
    do i_b = 1,b_num
        do i_y = 1,y_num
            READ(10, '(F15.10, X, F15.10, X, F15.10, X, F15.10)') v_matrix(i_b, i_y), &
                v0_matrix(i_b, i_y), v1_matrix(i_b, i_y), v_excl_matrix(i_b, i_y)
            READ(16, '(F15.11, X, F15.11)') q_matrix(i_b, i_y), q_matrix_nodef(i_b, i_y)
            READ(13, '(F15.11)') b_next_matrix(i_b, i_y)
            if (v1_matrix(i_b, i_y) > (v0_matrix(i_b, i_y))) then
                default_decision(i_b, i_y) =2
            else
                default_decision(i_b, i_y) =1
            end if
        end do
    end do
CLOSE(10)
!CLOSE(12)
CLOSE(16)
CLOSE(13)
CLOSE(17)
end if

call iterate
!call simulate
!call compute_initial_guess

call cpu_time(end_time)
WRITE(nout, '(A7, X, A7, X, A7)') 'Hours ', 'Minutes', 'Seconds'
WRITE(nout, '(I7, X, I7, X, I7)') INT((end_time - start_time) / 3600d+0), &
             INT((end_time-start_time)/60d+0 - INT((end_time - start_time) / 3600d+0)*60d+0),&
INT(end_time-start_time - INT((end_time - start_time) / 3600d+0)*3600d+0 - &
INT((end_time-start_time)/60d+0 - INT((end_time - start_time) / 3600d+0)*60d+0)*60d+0)
end program
