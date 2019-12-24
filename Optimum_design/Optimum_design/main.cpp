#include <stdio.h>
#include <math.h>
#include <Windows.h>
#include <stdlib.h>
#include <string.h>

#define DIM 3							// Dimension array 선언에 들어갈것들
#define INEQUALITY 6					// Equality Constraints
#define EQUALITY 2						// Inequality Constraints
//#define PHI sqrt(2)						// φ = 2^(1/3) case
//#define PHI pow(2,1.0/3.0) 
#define GOLDENRATIO (sqrt(5) + 1) / 2	// 1.618

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))  
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))  

// 식별자는 전체 대문자
// 함수이름은 첫글자 대문자
// 변수이름은 모두 소문자

long double PHI;
void Quasi_Newton_method(long double x[], int m);
long double Objective_function(long double x[]);
long double Penalty_function(long double x[], long double r);
void Gradient(long double *g, long double x0[], long double r);
void Hessian_BFGS(long double(*h)[DIM], long double a, long double d[], long double g1[], long double g0[]);
long double Norm(long double *x, int n); //Norm(vector, sizeof(vector) / sizeof(long double))
void Identity_Matrix(long double(*h)[DIM]);
long double Golden_section_method(long double X[], long double d[], long double r);
long double Quadraticapprox_withGolden(long double X[], long double d[], long double r);
long double Equalintervalsection(long double X[], long double d[], long double r);
void Inverse3(long double(*xi)[DIM], long double x[][DIM]);
void LU_decomposition(long double *d, long double h[][DIM], long double g[]);

int allstep = 0;
int rstep = 0; int sumstep = 0;

int main() {

	//system("color f0");
	long double design_variable[DIM];
	int onedmode;
	//printf("  1=> Golden Section Method,  2=> Quadratic Method  ?\n");
	//scanf_s("%d", &onedmode);
	fflush(stdin);
	int i_phi;
	
	for (int i = 0; i < 2; i++) {
		if (i == 0) {
			printf("  PHI --> 2^(1/2) \n\n");
			PHI = pow(2, 1.0 / 2.0);
		}
		else {
			printf("  PHI --> 2^(1/3) \n\n");
			PHI = pow(2, 1.0 / 3.0);
		}

		onedmode = 2; // 1 : Golden section,  2 : Quadratic approx. with golden section
		for (onedmode = 1; onedmode <= 2; onedmode++) {
			design_variable[0] = 21; design_variable[1] = 21; design_variable[2] = 21; // Initialize design variables
			Quasi_Newton_method(design_variable, onedmode);
		}

	}
	getchar();
	getchar();
	return 0;
}

long double Objective_function(long double x[]) {
	// Design variable 값을 받아서 Objective function 값을 반환
	return (x[1] - x[2])*(x[1] - x[2]);
}

long double Penalty_function(long double x[], long double r) {
	long double h[EQUALITY];
	long double g[INEQUALITY];
	long double p = 0;

	//h[0] = pow(PHI, 2)*x[0] * (x[0] + x[2] - x[1]) - x[1] * x[2];
	//h[1] = pow(PHI, 3)*x[0] - x[1] * (1 + x[1] - x[0]);
	h[0] = pow(PHI, 2) * x[0] * (x[0] + x[2] - x[1]) / (x[1] * x[2]) - 1;
	h[1] = 1 - x[1] * (1 + x[1] - x[0]) / (pow(PHI, 3)*x[0]);
	g[0] = x[0] / 20 - 1;
	g[1] = x[1] / 20 - 1;
	g[2] = x[2] / 20 - 1;
	g[3] = -x[0];
	g[4] = -x[1];
	g[5] = -x[2];

	for (int i = 0; i < EQUALITY; i++) {
		p += h[i] * h[i];
	}
	for (int i = 0; i < INEQUALITY; i++) {
		g[i] = MAX(g[i], 0);
		p += g[i] * g[i];
	}

	return Objective_function(x) + r * p;
}

void Gradient(long double *g, long double x0[], long double r) {

	// Design variable 값을 받아서 Gradient를 반환
	long double dx = 1e-7;
	long double fx = Penalty_function(x0, r);
	long double x1[DIM];

	// memcpy(x1, x0, sizeof(x0));
	for (int i = 0; i < DIM; i++) {
		x1[i] = x0[i];
	}
	
	for (int i = 0; i < DIM; i++) {
		x1[i] += dx;
		g[i] = (Penalty_function(x1, r) - fx) / dx;
		x1[i] = x0[i];
	}
}
//(hessian, design_variable, step_size, direction, gradient, gradient_temp)
void Hessian_BFGS(long double(*h)[DIM], long double a, long double d[], long double g1[], long double g0[]) {
	long double y[DIM];
	long double s[DIM];
	long double D[DIM][DIM];
	long double E[DIM][DIM];
	long double ydots = 0;
	long double cdotd = 0;

	for (int i = 0; i < DIM; i++) {
		y[i] = g1[i] - g0[i];			// change in gradient
		s[i] = a * d[i];				// change in design
		ydots += y[i] * s[i];
		cdotd += g0[i] * d[i];
	}

	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			D[i][j] = (y[i] * y[j]) / ydots;
			E[i][j] = (g0[i] * g0[j]) / cdotd;
			h[i][j] += D[i][j] + E[i][j];
			//printf("%lf , %lf , %lf\n", ydots, D[i][j], E[i][j]);
		}
	}
}

long double Norm(long double *x, int n) {
	long double accum = 0;
	for (int i = 0; i < n; i++) {
		accum += x[i] * x[i];
	}
	return sqrt(accum);
}

void Identity_Matrix(long double(*h)[DIM]) {
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			if (i == j) h[i][j] = 1;
			else h[i][j] = 0;
		}
	}
}

void Quasi_Newton_method(long double design_variable[], int m) {
	
	
	long double change_in_design[DIM];
	long double step_size;
	long double direction[DIM] = { 0, };
	long double penalty_parameter;
	long double termination_criterion;
	long double hessian[DIM][DIM];
	long double gradient[DIM];
	long double gradient_temp[DIM];
	long double hessian_inverse[DIM][DIM];
	long double norm;
	char onedmode_name[40];
	int iter = 1;
	int penalty_parameter_flag;

	if (m == 1) strcpy_s(onedmode_name, sizeof("Golden section method"), "Golden section method");
	else if (m ==2) strcpy_s(onedmode_name, sizeof("Golden section with Quadratic approx."), "Golden section with Quadratic approx.");
	else strcpy_s(onedmode_name, sizeof("Equal-interval section method"), "Equal-interval section method");

	Identity_Matrix(hessian);				// DIM 3, Identity, Initialize		
	termination_criterion = 1e-5;			// Norm 이 일정 크기 이하가 되면 Loop를 중단한다. 
	penalty_parameter = 1e-4;				
	//penalty_parameter = 10;
	
	printf("|  Quasi-Newton Method   %-40s            Termination criterion = %.0e            |\n", onedmode_name, termination_criterion);
	printf("---ITER---STEP SIZE-------NORM----------------------------DESIGN VARIABLES---------------------------------------------\n");
	for (int k = 0; k < 6; k++) {				// 초기 penalty parameter는 1e-5에서 시작하여 10배씩 증가시킴, 현재는 10으로 고정되어있음 (10으로 고정되어 현재 내부 Loop인 while문을 한번 나오게 되면 더이상 for loop는 돌지 않음) 
		penalty_parameter_flag = 1;
		Gradient(gradient, design_variable, penalty_parameter);
		Identity_Matrix(hessian);				// DIM 3, Identity, Initialize	
		norm = 1.;
		iter = 1;
		//norm = Norm(gradient, DIM);
		while (norm > termination_criterion) {
			
			// termination_criterion 을 penalty_parameter 에 대해 적절한 크기로 정하도록 만들기

			// Hessian inverse를 구하여 -gradient와 곱하는 방법
			// LU decomposition 이용하는 방법
			Inverse3(hessian_inverse, hessian);
			for (int i = 0; i < DIM; i++) {
				direction[i] = 0; //이전 direction 값 초기화
				for (int j = 0; j < DIM; j++) {
					direction[i] += hessian_inverse[i][j] * -gradient[j];
				}
			}
			//LU_decomposition(direction, hessian, gradient);
			
			// 결과 출력부분 Penalty parameter가 커질때마다 새로 작성
			if (penalty_parameter_flag) {
				printf("------------------------------------------- Penalty parameter = %1.2e ----------------------------------------------\n", penalty_parameter);
				penalty_parameter_flag = 0;
			}
			if(m == 1) step_size = Golden_section_method(design_variable, direction, penalty_parameter);		// Golden section method로 step size 정하기
			else if (m==2)step_size = Quadraticapprox_withGolden(design_variable, direction, penalty_parameter);				// Quadraticapproximation method로 step size 정하기
			else step_size =Equalintervalsection(design_variable, direction, penalty_parameter);

			for (int i = 0; i < DIM; i++) {
				gradient_temp[i] = gradient[i];														// Old gradient를 임시로 저장, BFGS에서 사용
				change_in_design[i] = step_size * direction[i];
				design_variable[i] += change_in_design[i];
			}
			// design_variable이 step size * direction으로 update됨
			Gradient(gradient, design_variable, penalty_parameter);									// Update한 Design variable에서 Gradient구하기
			Hessian_BFGS(hessian, step_size, direction, gradient, gradient_temp);					// hessian을 BFGS방법으로 Update, 필요한것 Gradient Old, New(Change in design, Change in gradient, Step size, Direction)
			norm = Norm(gradient, DIM);																// hessian_BFGS에서 -nan(ind) 0을 나누는듯

			//printf("%10.5e  %10.5e  %10.5e  %10.5e\n", step_size, direction, gradient, gradient_temp);
			printf("|  %4d  %10.5e  %10.5e  %20.17lf  %20.17lf  %20.17lf                   |\n", iter, step_size, norm, design_variable[0], design_variable[1], design_variable[2]);

			iter++;
			// getchar();
		}
		//termination_criterion *= 0.1;
		penalty_parameter *= 10;
		
		allstep;
		printf("\nsteps for specific r : %d \n\n", allstep);
		sumstep += allstep;
		allstep = 0;
	}
	printf("|--------------------------------------------------------END----------------------------------------------------------|\n");
	printf("                         Sum of step = %d                         \n\n\n\n",sumstep);
	sumstep = 0; //초기화
}

void Inverse3(long double(*xi)[DIM], long double x[][DIM]) {
	long double det = 0;
	for (int i = 0; i < DIM; i++) det += (x[0][i] * (x[1][(i + 1) % 3] * x[2][(i + 2) % 3] - x[1][(i + 2) % 3] * x[2][(i + 1) % 3]));
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			xi[i][j] = ((x[(i + 1) % 3][(j + 1) % 3] * x[(i + 2) % 3][(j + 2) % 3]) - (x[(i + 1) % 3][(j + 2) % 3] * x[(i + 2) % 3][(j + 1) % 3])) / det;
		}
	}
}

long double Golden_section_method(long double X[], long double d[], long double r) {
	long double x0[DIM];
	long double x1[DIM];
	long double x2[DIM];
	long double alpha;
	long double interval;
	long double delta_ini = 0.5;
	long double epsilon = 0.0001;
	long double lower_bound, medium_lower, medium_upper, upper_bound;
	long double initial_parameter[DIM], y0, y1, y2;
	long double x_a[DIM], x_b[DIM];
	long double golden_ratio = GOLDENRATIO;
	int iteration = 0, i;
	int iteration_initial = 0;
	int step_number = 0;


	do {
		delta_ini *= 0.1;

		for (i = 0; i < DIM; i++) {
			initial_parameter[i] = (delta_ini * ((pow(golden_ratio, i + 1) - 1)) / (golden_ratio - 1));
		}

		for (i = 0; i < DIM; i++) {
			x0[i] = X[i] + initial_parameter[0] * d[i];
		}
		y0 = Penalty_function(x0, r);

		for (i = 0; i < DIM; i++) {
			x1[i] = X[i] + initial_parameter[1] * d[i];
		}
		y1 = Penalty_function(x1, r);
		step_number += 1;//newly step

	} while (y1 > y0);
	//printf("\n initial stepsize decision step numbers: %d", step_number);
	for (i = 0; i < DIM; i++) {
		x2[i] = X[i] + initial_parameter[2] * d[i];
	}
	y2 = Penalty_function(x2, r);
		
	//printf("%f %f %f\n", initial_parameter[0], initial_parameter[1], initial_parameter[2]);
	//printf("%f %f %f\n", x0[0], x0[1], x0[2]);
	//printf("%f %f %f\n", x1[0], x1[1], x1[2]);
	//printf("%f %f %f\n", x2[0], x2[1], x2[2]);
	//printf("%.10lf %.10lf %.10lf\n", y0, y1, y2);
	//intializtion

	while (y1 > y2 || y1 > y0) {
		initial_parameter[0] = delta_ini * (pow(golden_ratio, iteration_initial + 2) - 1) / (golden_ratio - 1);
		initial_parameter[1] = delta_ini * (pow(golden_ratio, iteration_initial + 3) - 1) / (golden_ratio - 1);
		initial_parameter[2] = delta_ini * (pow(golden_ratio, iteration_initial + 4) - 1) / (golden_ratio - 1);

		for (i = 0; i < DIM; i++) {
			x0[i] = X[i] + initial_parameter[0] * d[i];
		}
		//y0 = penalty_function(x0, r);
		y0 = y1;


		for (i = 0; i < DIM; i++) {
			x1[i] = X[i] + initial_parameter[1] * d[i];
		}
		//y1 = penalty_function(x1, r);
		y1 = y2;

		for (i = 0; i < DIM; i++) {
			x2[i] = X[i] + initial_parameter[2] * d[i];
		}
		y2 = Penalty_function(x2, r);

		//printf("\ny0:%e y1:%e y2:%e\n",y0, y1, y2);
		iteration_initial++;
		step_number += 1;//newly step
	}
	//printf("\n                                                           Interval searching steps : %d , ", step_number);
	lower_bound = initial_parameter[0];
	upper_bound = initial_parameter[2];
	interval = upper_bound - lower_bound;
	medium_lower = lower_bound + (2 - golden_ratio)*interval;//already stepped 
	medium_upper = lower_bound + (golden_ratio - 1)*interval;//newly step.
	step_number += 1;//newly step

	if (lower_bound < 0) {
		goto back;
	}

	long double y_l, y_a, y_b, y_u;

	long double X_l[DIM];
	for (i = 0; i < DIM; i++) {
		X_l[i] = X[i] + lower_bound * d[i];
	}
	//y_l = penalty_function(X_l, r);
	y_l = y0;

	long double X_a[DIM];
	for (i = 0; i < DIM; i++) {
		X_a[i] = X[i] + medium_lower * d[i];
	}
	//y_a = penalty_function(X_a, r);
	y_a = y1;
	
	long double X_b[DIM];
	for (i = 0; i < DIM; i++) {
		X_b[i] = X[i] + medium_upper * d[i];
	}
	y_b = Penalty_function(X_b, r);

	long double X_u[DIM];
	for (i = 0; i < DIM; i++) {
		X_u[i] = X[i] + upper_bound * d[i];
	}
	//y_u = Penalty_function(X_u, r);
	y_u = y2;

	while (interval > epsilon) {

		if (y_a < y_b) {
			upper_bound = medium_upper;
			y_u = y_b;
			medium_upper = medium_lower;
			y_b = y_a;
			interval = upper_bound - lower_bound;
			medium_lower = lower_bound + (2 - golden_ratio) * interval;
			//X_a = Calculate_step(X, medium_lower, d);
			for (i = 0; i < DIM; i++) {
				X_a[i] = X[i] + medium_lower * d[i];
			}
			y_a = Penalty_function(X_a, r);
		}

		else {
			lower_bound = medium_lower;
			y_l = y_a;
			medium_lower = medium_upper;
			y_a = y_b;
			interval = upper_bound - lower_bound;
			medium_upper = lower_bound + (golden_ratio - 1) * interval;
			//X_b = Calculate_step(X, medium_upper, d);
			for (i = 0; i < DIM; i++) {
				x_b[i] = X[i] + medium_upper * d[i];
			}
			y_b = Penalty_function(x_b, r);
		}

		//printf("Golden_iteration: %d lower: %.10f a: %.10f b: %.10f upper: %.10f \n", iteration, lower_bound, medium_lower, medium_upper, upper_bound);

		iteration++;
		step_number += 1;//newly step
	}

	//X_a = Calculate_step(X, medium_lower, d);
	//X_b = Calculate_step(X, medium_upper, d);

	for (i = 0; i < DIM; i++) {
		X_a[i] = X[i] + medium_lower * d[i];
	}

	for (i = 0; i < DIM; i++) {
		X_b[i] = X[i] + medium_upper * d[i];		
	}

	alpha = (medium_lower + medium_upper) / 2;
	step_number += 1;//newly step
	//printf("Total Step numbers : %d                        \n", step_number);
	allstep += step_number;
	//printf("\n\n alpha : %lf", alpha);
back:
	return alpha;
	
}

long double Equalintervalsection(long double X[], long double d[], long double r) {

	long double x0[DIM];
	long double x1[DIM];
	long double x2[DIM];
	long double alpha;
	long double interval;
	long double delta_ini = 0.5;
	long double epsilon = 0.0001;
	long double lower_bound, medium_lower, medium_upper, upper_bound;
	long double initial_parameter[DIM], y0, y1, y2;
	long double x_a[DIM], x_b[DIM];
	int iteration = 0, i;
	int iteration_initial = 0;
	int step_number = 0;
	long double e1= 1.0 / 3.0;
	long double e2 = 2.0 / 3.0;


	do {
		delta_ini *= 0.1;

		for (i = 0; i < DIM; i++) {
			initial_parameter[i] = delta_ini * (i);
		}
		//printf("%f %f %f", initial_parameter[0], initial_parameter[1], initial_parameter[2]);


		for (i = 0; i < DIM; i++) {
			x0[i] = X[i] + initial_parameter[0] * d[i];
		}
		y0 = Penalty_function(x0, r);


		for (i = 0; i < DIM; i++) {
			x1[i] = X[i] + initial_parameter[1] * d[i];
		}
		y1 = Penalty_function(x1, r);
		step_number += 1;//newly step
	} while (y1 > y0);
	//printf("\n initial stepsize decision step numbers: %d", step_number);
	
	for (i = 0; i < DIM; i++) {
		x2[i] = X[i] + initial_parameter[2] * d[i];
	}
	y2 = Penalty_function(x2, r);


	//printf("%f %f %f\n", initial_parameter[0], initial_parameter[1], initial_parameter[2]);
	//printf("%f %f %f\n", x0[0], x0[1], x0[2]);
	//printf("%f %f %f\n", x1[0], x1[1], x1[2]);
	//printf("%f %f %f\n", x2[0], x2[1], x2[2]);
	//printf("%.10lf %.10lf %.10lf\n", y0, y1, y2);
	//intializtion
	while (y1 > y2 || y1 > y0) {
		initial_parameter[0] = delta_ini * (iteration_initial + 1);
		initial_parameter[1] = delta_ini * (iteration_initial + 2);
		initial_parameter[2] = delta_ini * (iteration_initial + 3);


		for (i = 0; i < DIM; i++) {
			x0[i] = X[i] + initial_parameter[0] * d[i];
		}
		//y0 = penalty_function(x0, r);
		y0 = y1;

		for (i = 0; i < DIM; i++) {
			x1[i] = X[i] + initial_parameter[1] * d[i];
		}
		//y1 = penalty_function(x1, r);
		y1 = y2;

		for (i = 0; i < DIM; i++) {
			x2[i] = X[i] + initial_parameter[2] * d[i];
		}
		y2 = Penalty_function(x2, r);

		//printf("\ny0:%f y1:%f y2:%f\n",y0, y1, y2);
		iteration_initial++;
		//printf("%.10lf %.10lf %.10lf\n", y0, y1, y2);
		step_number += 1;//newly step
	}
	//printf("\n                                                           Interval searching steps : %d , ", step_number);
	//printf("Equal_iteration_initial: %d \n", iteration_initial);

	lower_bound = initial_parameter[0];
	upper_bound = initial_parameter[2];
	interval = upper_bound - lower_bound;
	medium_lower = lower_bound + (e1)*interval;
	medium_upper = lower_bound + (e2)*interval;
	if (lower_bound < 0) {
		goto back;
	}

	iteration = 0;

	long double y_l, y_a, y_b, y_u;

	//long double * X_l = Calculate_step(X, lower_bound, d);

	long double X_l[DIM];
	for (i = 0; i < DIM; i++) {
		X_l[i] = X[i] + lower_bound * d[i];
	}

	//y_l = penalty_function(X_l, r);
	y_l = y0;

	//long double * X_a = Calculate_step(X, medium_lower, d);
	long double X_a[DIM];
	for (i = 0; i < DIM; i++) {
		X_a[i] = X[i] + medium_lower * d[i];
	}
	y_a = Penalty_function(X_a, r);
	//y_a = y1;
	//long double * X_b = Calculate_step(X, medium_upper, d);

	long double X_b[DIM];
	for (i = 0; i < DIM; i++) {
		X_b[i] = X[i] + medium_upper * d[i];
	}
	y_b = Penalty_function(X_b, r);

	//long double * X_u = Calculate_step(X, upper_bound, d);
	long double X_u[DIM];
	for (i = 0; i < DIM; i++) {
		X_u[i] = X[i] + medium_upper * d[i];
	}
	//y_u = penalty_function(X_u, r);
	y_u = y2;

	while (interval > epsilon) {

		if (y_a < y_b) {
			upper_bound = medium_upper;
			y_u = y_b;
			interval = upper_bound - lower_bound;
		}

		else {
			lower_bound = medium_lower;
			y_l = y_a;
			interval = upper_bound - lower_bound;
		}

		medium_lower = lower_bound + (e1)* interval;
		//X_a = Calculate_step(X, medium_lower, d);
		for (i = 0; i < DIM; i++) {
			X_a[i] = X[i] + medium_lower * d[i];
		}
		y_a = Penalty_function(X_a, r);

		medium_upper = lower_bound + (e2)* interval;
		for (i = 0; i < DIM; i++) {
			X_b[i] = X[i] + medium_upper * d[i];
		}
		y_b = Penalty_function(X_b, r);;
		step_number += 2;//newly 2 steps
		
		iteration +=2;//equal-interval 에서는 1 iteration 마다 1 points가 아닌 2 points가 찍히므로 횟수에 고려해야 맞음.
		//printf("Equal_iteration: %d lower: %.10f a: %.10f b: %.10f upper: %.10f \n", iteration, lower_bound, medium_lower, medium_upper, upper_bound);

	}

	//X_a = Calculate_step(X, medium_lower, d);
	//X_b = Calculate_step(X, medium_upper, d);

	for (i = 0; i < DIM; i++) {
		X_a[i] = X[i] + medium_lower * d[i];
	}

	for (i = 0; i < DIM; i++) {
		X_b[i] = X[i] + medium_upper * d[i];
	}

	alpha = (medium_lower + medium_upper) / 2;
	step_number += 1;//newly step
	//printf("Total Step numbers : %d                        \n", step_number);
	allstep += step_number;
	//printf("\n\n alpha : %lf", alpha);
back:
	return alpha;
}

long double Quadraticapprox_withGolden(long double X[], long double d[], long double r) {

	long double x0[DIM];
	long double x1[DIM];
	long double x2[DIM];
	long double alpha;
	long double interval;
	long double delta_ini = 0.5;
	long double epsilon = 0.0001;
	long double lower_bound, medium_lower, upper_bound, medium_i, medium_bar, medium_bar_old;
	long double initial_parameter[DIM], y0, y1, y2;
	long double dybar = 1, dxbar = 100;
	long double golden_ratio = GOLDENRATIO;
	int iteration = 0, i;
	int iteration_initial = 0;
	int step_number = 0;

	do {
		delta_ini *= 0.1;

		for (i = 0; i < DIM; i++) {
			initial_parameter[i] = (delta_ini * ((pow(golden_ratio, i + 1) - 1)) / (golden_ratio - 1));
		}

		for (i = 0; i < DIM; i++) {
			x0[i] = X[i] + initial_parameter[0] * d[i];
		}
		y0 = Penalty_function(x0, r);

		for (i = 0; i < DIM; i++) {
			x1[i] = X[i] + initial_parameter[1] * d[i];
		}
		y1 = Penalty_function(x1, r);
		step_number += 1;//newly step
	} while (y1 > y0);
	//printf("\n initial stepsize decision step numbers: %d", step_number);

	for (i = 0; i < DIM; i++) {
		x2[i] = X[i] + initial_parameter[2] * d[i];
	}
	y2 = Penalty_function(x2, r);
	
	//printf("initial: %f %f %f\n", initial_parameter[0], initial_parameter[1], initial_parameter[2]);
	//printf("%f %f %f\n", x0[0], x0[1], x0[2]);
	//printf("%f %f %f\n", x1[0], x1[1], x1[2]);
	//printf("%f %f %f\n", x2[0], x2[1], x2[2]);
	//printf("%.10lf %.10lf %.10lf\n", y0, y1, y2);
	//intializtion
	while (y1 > y2 || y1 > y0) {
		initial_parameter[0] = delta_ini * (pow(golden_ratio, iteration_initial + 2) - 1) / (golden_ratio - 1);
		initial_parameter[1] = delta_ini * (pow(golden_ratio, iteration_initial + 3) - 1) / (golden_ratio - 1);
		initial_parameter[2] = delta_ini * (pow(golden_ratio, iteration_initial + 4) - 1) / (golden_ratio - 1);

		for (i = 0; i < DIM; i++) {
			x0[i] = X[i] + initial_parameter[0] * d[i];
		}
		//y0 = penalty_function(x0, r);
		y0 = y1;


		for (i = 0; i < DIM; i++) {
			x1[i] = X[i] + initial_parameter[1] * d[i];
		}
		//y1 = penalty_function(x1, r);
		y1 = y2;

		for (i = 0; i < DIM; i++) {
			x2[i] = X[i] + initial_parameter[2] * d[i];
		}
		y2 = Penalty_function(x2, r);
		step_number += 1;//newly step
		//printf("\ny0:%f y1:%f y2:%f\n",y0, y1, y2);
		iteration_initial++;
	}
	//printf("\n                                                           Interval searching steps : %d , ", step_number);
	lower_bound = initial_parameter[0];
	upper_bound = initial_parameter[2];
	interval = upper_bound - lower_bound;
	medium_i = lower_bound + (2 - golden_ratio) * interval;

	long double y_bar_old = 0;
	long double y_l, y_i, y_bar, y_u;
	if (lower_bound < 0) {
		goto back;
	}

	iteration = 0;
	
	long double X_l[DIM];
	for (i = 0; i < DIM; i++) {
		X_l[i] = X[i] + lower_bound * d[i];
	}

	//y_l = penalty_function(X_l, r);
	y_l = y0;

	long double X_i[DIM];
	for (i = 0; i < DIM; i++) {
		X_i[i] = X[i] + medium_i * d[i];
	}
	//y_i = penalty_function(X_i, r);
	y_i = y1;

	/*
	long double X_b[Dim];
	for (i = 0; i < Dim; i++) {
		X_b[i] = X[i] + medium_upper * d[i];
	}
	y_b = penalty_function(X_b, r);
	*/

	long double X_u[DIM];
	for (i = 0; i < DIM; i++) {
		X_u[i] = X[i] + upper_bound * d[i];
	}
	//y_u = penalty_function(X_u, r);
	y_u = y2;

	long double a0, a1, a2;
	long double a0s, a1s, a2s, medium_bars;
	a2 = (1 / (upper_bound - medium_i))*((y_u - y_l) / (upper_bound - lower_bound) - (y_i - y_l) / (medium_i - lower_bound));
	a1 = (y_i - y_l) / (medium_i - lower_bound) - a2 * (lower_bound + medium_i);
	a0 = y_l - (a1*lower_bound) - a2 * pow(lower_bound, 2.0);
	medium_bar = -(a1) / (2.0 * a2);//necessary condition: a2>0 which means "u" form.
	step_number += 1;//newly step

	long double X_bar[DIM];
	for (i = 0; i < DIM; i++) {
		X_bar[i] = X[i] + medium_bar * d[i];
	}
	y_bar = Penalty_function(X_bar, r);

	a2s = a2; a1s = a1; a0s = a0;
	medium_bars = medium_bar;
	//printf("\n before while a2: %10.5e , a1: %10.5e , a0: %10.5e , medium_bar %10.5e \n", a2, a1, a0, medium_bar);
	while (interval > epsilon && a2>0.0 && a1<0.0) {
		if (medium_i < medium_bar) {
			if (y_i < y_bar) {
				//printf("\n 1 \n");
				upper_bound = medium_bar;
				y_u = y_bar;
				y_bar_old = y_bar;
				medium_bar_old = medium_bar;

				a2 = (1 / (upper_bound - medium_i))*((y_u - y_l) / (upper_bound - lower_bound) - (y_i - y_l) / (medium_i - lower_bound));
				a1 = (y_i - y_l) / (medium_i - lower_bound) - a2 * (lower_bound + medium_i);
				a0 = y_l - (a1*lower_bound) - a2 * pow(lower_bound, 2);
				medium_bar = -(a1) / (2 * a2);//necessary condition: a2>0 which means "u" form.
				

				long double X_bar[DIM];
				for (i = 0; i < DIM; i++) {
					X_bar[i] = X[i] + medium_bar * d[i];
				}
				y_bar = Penalty_function(X_bar, r);

				dybar = (y_bar - y_bar_old);
				dxbar = (medium_bar - medium_bar_old);
				if (dybar < 0) {
					dybar = -dybar;
				}
				if (dxbar < 0) {
					dxbar = -dxbar;
				}
				interval = upper_bound - lower_bound;
			}

			else {
				//printf("\n 2 \n");
				lower_bound = medium_i;
				y_l = y_i;
				medium_i = medium_bar;
				y_i = y_bar;
				medium_bar_old = medium_bar;

				a2 = (1 / (upper_bound - medium_i))*((y_u - y_l) / (upper_bound - lower_bound) - (y_i - y_l) / (medium_i - lower_bound));
				a1 = (y_i - y_l) / (medium_i - lower_bound) - a2 * (lower_bound + medium_i);
				a0 = y_l - (a1*lower_bound) - a2 * pow(lower_bound, 2);
				medium_bar = -(a1) / (2 * a2);//necessary condition: a2>0 which means "u" form.

				long double X_bar[DIM];
				for (i = 0; i < DIM; i++) {
					X_bar[i] = X[i] + medium_bar * d[i];
				}
				y_bar = Penalty_function(X_bar, r);

				dybar = (y_bar - y_bar_old);
				dxbar = (medium_bar - medium_bar_old);
				if (dybar < 0) {
					dybar = -dybar;
				}
				if (dxbar < 0) {
					dxbar = -dxbar;
				}
				interval = upper_bound - lower_bound;
			}
		}
		else {
			if (y_i < y_bar) {
				//printf("\n 3 \n");
				lower_bound = medium_bar;
				y_l = y_bar;
				y_bar_old = y_bar;
				medium_bar_old = medium_bar;

				a2 = (1 / (upper_bound - medium_i))*((y_u - y_l) / (upper_bound - lower_bound) - (y_i - y_l) / (medium_i - lower_bound));
				a1 = (y_i - y_l) / (medium_i - lower_bound) - a2 * (lower_bound + medium_i);
				a0 = y_l - (a1*lower_bound) - a2 * pow(lower_bound, 2);
				medium_bar = -(a1) / (2 * a2);//necessary condition: a2>0 which means "u" form.

				long double X_bar[DIM];
				for (i = 0; i < DIM; i++) {
					X_bar[i] = X[i] + medium_bar * d[i];
				}
				y_bar = Penalty_function(X_bar, r);

				dybar = (y_bar - y_bar_old);
				dxbar = (medium_bar - medium_bar_old);
				if (dybar < 0) {
					dybar = -dybar;
				}
				if (dxbar < 0) {
					dxbar = -dxbar;
				}
				interval = upper_bound - lower_bound;
			}

			else {
				//printf("\n 4 \n");
				upper_bound = medium_i;
				y_u = y_i;
				medium_i = medium_bar;
				y_i = y_bar;
				medium_bar_old = medium_bar;

				a2 = (1 / (upper_bound - medium_i))*((y_u - y_l) / (upper_bound - lower_bound) - (y_i - y_l) / (medium_i - lower_bound));
				a1 = (y_i - y_l) / (medium_i - lower_bound) - a2 * (lower_bound + medium_i);
				a0 = y_l - (a1*lower_bound) - a2 * pow(lower_bound, 2);
				medium_bar = -(a1) / (2 * a2);//necessary condition: a2>0 which means "u" form.

				long double X_bar[DIM];
				for (i = 0; i < DIM; i++) {
					X_bar[i] = X[i] + medium_bar * d[i];
				}
				y_bar = Penalty_function(X_bar, r);

				dybar = (y_bar - y_bar_old);
				dxbar = (medium_bar - medium_bar_old);
				if (dybar < 0){
					dybar = -dybar;
				}
				if (dxbar < 0) {
					dxbar = -dxbar;
				}
				interval = upper_bound - lower_bound;
			}
			

		}

		if (a2 > 0.0) {
			a2s = a2; a1s = a1; a0s = a0; medium_bars = medium_bar;
			iteration++;
			step_number += 1;//newly step
			//printf("\n a2: %10.5e , a1: %10.5e , a0: %10.5e , medium_bar %10.5e \n", a2s, a1s, a0s, medium_bars);
		}
		
		//printf("quadratic_iteration: %d lower: %.10f ui: %.10f ubar: %.10f upper: %.10f \n", iteration, lower_bound, medium_i, medium_bar, upper_bound);
		//printf("\n a2: %f, a1: %f, a0= %f\n", a2, a1, a0);
		
	}

	//X_a = Calculate_step(X, medium_lower, d);
	//X_b = Calculate_step(X, medium_upper, d);

	alpha = medium_bars;
	step_number += 1;//newly step
	//printf("Total Step numbers : %d                        \n", step_number);
	allstep += step_number;
	//printf("\n\n alpha : %lf", alpha);
back:
	return alpha;
}

void LU_decomposition(long double *d, long double h[][DIM], long double g[]) {
	long double temp_matrix[DIM][DIM];
	long double L[DIM][DIM];
	long double y[DIM];
	long double sum = 0;
	long double pivot;
	static long double X[DIM];
	int i = 0, j = 0, k = 0;

	memcpy(temp_matrix, h, sizeof(long double)*DIM*DIM);

	for (i = 0; i < DIM; i++) {
		for (j = 0; j < DIM; j++) {
			if (i == j) L[i][j] = 1;
			else L[i][j] = 0;
		}
	}

	for (i = 1; i < DIM; i++) {
		if (temp_matrix[i - 1][i - 1] != 0) {
			for (j = i; j < DIM; j++) {
				if (temp_matrix[j][i] == 0) pivot = 0;
				else {
					pivot = temp_matrix[j][i - 1] / temp_matrix[i - 1][i - 1];
					for (k = 0; k < DIM; k++) temp_matrix[j][k] -= pivot * temp_matrix[i - 1][k];
				}
				L[j][i - 1] = pivot;
			}
		}
	}

	for (i = 0; i < DIM; i++) {
		for (j = 0; j < i; j++) {
			sum -= L[i][j] * y[j];
		}
		y[i] = g[i] + sum;
		sum = 0;
	}

	sum = 0;
	for (i = DIM - 1; i >= 0; i--) {
		for (j = DIM - 1; j > i; j--) {
			sum -= temp_matrix[i][j] * X[j];
		}
		d[i] = (y[i] + sum) / temp_matrix[i][i];
		sum = 0;
	}
}
// printf("%e \n", x0[i]); //Sleep(500);

