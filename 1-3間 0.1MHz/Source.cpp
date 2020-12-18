#include<stdio.h>
#include<math.h>
#include<time.h>

//ループ回数(x,y,z)//
int Nx = 40;
int Ny = 251;
int Nz = 40;

/*解析領域*/
//x=120[mm],y=753[mm],z=120[mm]

//x方向の導線
#define x_s  20			//x方向導線
//y方向の導線
#define y_s  17			//y方向導線始まり(左端から51[mm]の場所)
#define y_5  164		//y方向(磁界測定場所　左端から51＋441[mm])
#define y_e  234		//y方向導線終わり(右端から124[mm]の場所)
//z方向の導線
#define z_s  17			//z方向導線始まり
#define z_v  20			//z方向導線～電源
#define z_e  23			//z方向導線終わり
#define z_r  17			//z方向導線～抵抗
#define z_re 23			//z方向抵抗～導線
//抵抗の位置
#define z_lr 21			//抵抗の場所
//インダクタンスの位置
#define z_ll 19			//インダクタンスの場所

//π//
double pi = 4.0 * atan(1.0);

//サイズ//
double dx = 3.0e-3;
double dy = 3.0e-3;
double dz = 3.0e-3;
double dt = 5.0e-12;
double ds = 3.0e-3;

//パラメータの値//
double c = 3.0e8;					//電磁波の進行速度　c[m/sec]
double σ = 1.0e-8;					//導電率　σ[S/m]
double μ = 4.0e-7*pi;				//透磁率　μ[H/m]
double ε = 1.0e-9 / (36 * pi);		//誘電率　ε[F/m]
double t;
double R = 1370;					//抵抗	R[Ω]

double L = 5.94e-3;					//インダクタンス L[H]

									//電界,磁界//
double Ex[41][252][41], Ex_old[41][252][41];
double Ey[41][252][41], Ey_old[41][252][41];
double Ez[41][252][41], Ez_old[41][252][41];

double Hx[41][252][41], Hx_old[41][252][41];
double Hy[41][252][41], Hy_old[41][252][41];
double Hz[41][252][41], Hz_old[41][252][41];

//電圧,電流
double V = 0.0;
double V1 = 0.0;
double Vf = 0.0;
double I = 0.0;
double H = 0.0;

//導線の設定
double r = 1.0e-3;			//導線の半径[m]
double r0 = 0.230;			//導線の等価半径[m]
double d = 1.8e-2;			//平行2線の距離[m]

							//等価媒質定数法
double m = 0;				//等価媒質定数
double εm = 0;				//導線近傍でのεm[F/m]
double μm = 0;				//導線近傍でのμm[H/m]


							//電界の計算の係数//
double E1 = 0;		double E2 = 0;
//等価媒質定数法を適当した時の電界の計算の係数//
double Em1 = 0;		double Em2 = 0;

//抵抗の電界の計算の係数//
double R1 = 0;      double R2 = 0;

//インダクタンスの電界の計算の係数//
double integral = 0;
//インダクタンスの理論値//
double L_t = 0;
double L_v = 0;
double L_i = 0;

int T1, T2, T3;
int i, j, k;
int count = 0;
char moji[101];

void main() {

	FILE *fp;
	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	FILE *fp4;
	FILE *fp5;

	sprintf_s(moji, "電圧の時間変化(im).csv", count);
	fopen_s(&fp, moji, "w");

	sprintf_s(moji, "電流の時間変化(im).csv", count);
	fopen_s(&fp1, moji, "w");

	sprintf_s(moji, "磁界観測点の場所の電圧(im).csv", count);
	fopen_s(&fp2, moji, "w");

	sprintf_s(moji, "磁界観測点の場所の電圧(im).txt", count);
	fopen_s(&fp3, moji, "w");

	sprintf_s(moji, "電流の時間変化(im).txt", count);
	fopen_s(&fp4, moji, "w");

	sprintf_s(moji, "磁界の時間変化(im).txt", count);
	fopen_s(&fp5, moji, "w");

	m = (log(1 / r0) / log(ds / r));		//等価媒質定数
	εm = ε*m;							    //導線近傍でのεm[F/m]
	μm = μ / m;						    //導線近傍でのμm[H/m]

	T1 = clock();

	for (i = 0; i <= Nx; i++) {
		for (j = 0; j <= Ny; j++) {
			for (k = 0; k <= Nz; k++) {
				Ex[i][j][k] = 0.0;	Ex_old[i][j][k] = 0.0;
				Ey[i][j][k] = 0.0;	Ey_old[i][j][k] = 0.0;
				Ez[i][j][k] = 0.0;	Ez_old[i][j][k] = 0.0;
				Hx[i][j][k] = 0.0;	Hx_old[i][j][k] = 0.0;
				Hy[i][j][k] = 0.0;	Hy_old[i][j][k] = 0.0;
				Hz[i][j][k] = 0.0;	Hz_old[i][j][k] = 0.0;
			}
		}
	}

	E1 = (2 * ε - σ*dt) / (2 * ε + σ*dt);			E2 = (2 * dt) / (σ*dt + 2 * ε);			//電界の計算の係数//
	Em1 = (2 * εm - σ*dt) / (2 * εm + σ*dt);		Em2 = (2 * dt) / (σ*dt + 2 * εm);			//等価媒質定数法を適当した時の電界の計算の係数//


	R1 = (2 * R*ε*dx*dy - dt * dz) / (2 * R*ε*dx*dy + dt * dz);
	R2 = dt*(2 * R*dx*dy) / (2 * ε*R*dx*dy + dz*dt);						//抵抗の電界の計算の係数//

	for (t = 0; t <= 1.0e-6; ) {

		//電界の計算//
		for (i = 0; i < Nx; i++) {
			for (j = 1; j < Ny; j++) {
				for (k = 1; k < Nz; k++) {

					Ex[i][j][k] = E1 * Ex_old[i][j][k] + (E2 / dy) * (Hz_old[i][j][k] - Hz_old[i][j - 1][k]) - (E2 / dz) * (Hy_old[i][j][k] - Hy_old[i][j][k - 1]);

				}
			}
		}
		for (i = 1; i < Nx; i++) {
			for (j = 0; j < Ny; j++) {
				for (k = 1; k < Nz; k++) {

					Ey[i][j][k] = E1 * Ey_old[i][j][k] + (E2 / dz) * (Hx_old[i][j][k] - Hx_old[i][j][k - 1]) - (E2 / dz) * (Hz_old[i][j][k] - Hz_old[i - 1][j][k]);

				}
			}
		}
		for (i = 1; i < Nx; i++) {
			for (j = 1; j < Ny; j++) {
				for (k = 0; k < Nz; k++) {

					Ez[i][j][k] = E1 * Ez_old[i][j][k] + (E2 / dx) * (Hy_old[i][j][k] - Hy_old[i - 1][j][k]) - (E2 / dy) * (Hx_old[i][j][k] - Hx_old[i][j - 1][k]);

				}
			}
		}

		//抵抗の電界計算//
		Ez[x_s][y_e][z_lr] = R1*Ez_old[x_s][y_e][z_lr] + (R2 / dx)*(Hy_old[x_s][y_e][z_lr] - Hy_old[x_s - 1][y_e][z_lr]) - (R2 / dy)*(Hx_old[x_s][y_e][z_lr] - Hx_old[x_s][y_e - 1][z_lr]);

		//インダクタンスの電界計算//
		Ez[x_s][y_e][z_ll] = 1.0*Ez_old[x_s][y_e][z_ll] + (dt / ε)*1.0*(Hy_old[x_s][y_e][z_ll] - Hy_old[x_s - 1][y_e][z_ll]) - (dt / ε)*1.0*(Hx_old[x_s][y_e][z_ll] - Hx_old[x_s][y_e - 1][z_ll]) - integral;

		//インダクタンスの計算の係数//
		integral += dz*dt*dt / (ε*L*dx*dy)*Ez[x_s][y_e][z_ll];

		//等価媒質定数法(電界)//
		for (j = y_s + 1; j < y_e; j++) {
			Ex[x_s - 1][j][z_e] = Em1 * Ex_old[x_s - 1][j][z_e] + (Em2 / dy) * (Hz_old[x_s - 1][j][z_e] - Hz_old[x_s - 1][j - 1][z_e]) - (Em2 / dz) * (Hy_old[x_s - 1][j][z_e] - Hy_old[x_s - 1][j][z_e - 1]);
			Ex[x_s][j][z_e] = Em1 * Ex_old[x_s][j][z_e] + (Em2 / dy) * (Hz_old[x_s][j][z_e] - Hz_old[x_s][j - 1][z_e]) - (Em2 / dz) * (Hy_old[x_s][j][z_e] - Hy_old[x_s][j][z_e - 1]);

			Ez[x_s][j][z_e - 1] = Em1 * Ez_old[x_s][j][z_e - 1] + (Em2 / dx) * (Hy_old[x_s][j][z_e - 1] - Hy_old[x_s - 1][j][z_e - 1]) - (Em2 / dy) * (Hx_old[x_s][j][z_e - 1] - Hx_old[x_s][j - 1][z_e - 1]);
			Ez[x_s][j][z_e] = Em1 * Ez_old[x_s][j][z_e] + (Em2 / dx) * (Hy_old[x_s][j][z_e] - Hy_old[x_s - 1][j][z_e]) - (Em2 / dy) * (Hx_old[x_s][j][z_e] - Hx_old[x_s][j - 1][z_e]);
		}

		//////////////////////////////////////////// 導線ｙ方向(上の導線) ////////////////////////////////////////////

		for (j = y_s; j < y_e; j++) {
			Ey[x_s][j][z_e] = 0.0;
		}

		//////////////////////////////////////////// 導線ｙ方向(下の導線) ////////////////////////////////////////////

		for (j = y_s; j < y_e; j++) {
			Ey[x_s][j][z_s] = 0.0;
		}

		//////////////////////////////////////////// 導線z方向 ////////////////////////////////////////////

		for (k = z_s; k < z_v; k++) {			//下の導線～電源
			Ez[x_s][y_s][k] = 0.0;
		}
		for (k = z_v + 1; k < z_e; k++) { 		//電源～上の導線
			Ez[x_s][y_s][k] = 0.0;
		}
		for (k = z_r; k < z_ll; k++) {			//下の導線～コイル
			Ez[x_s][y_e][k] = 0.0;
		}
		for (k = z_ll + 1; k < z_lr; k++) {		//コイル～抵抗
			Ez[x_s][y_e][k] = 0.0;
		}
		for (k = z_lr + 1; k < z_re; k++) {      //コイル～導線
			Ez[x_s][y_e][k] = 0.0;
		}

		//電界から電圧を与える式
		if (t <= 4.0e-8) {
			Ez[x_s][y_s][z_v] = (t*(0.2 / 4e-8)) / dz;
			V = Ez[x_s][y_s][z_v] * dz;
		}
		else if (t > 4.0e-8&&t <= 6.0e-8) {
			Ez[x_s][y_s][z_v] = (t*(0.64 / 2e-8) - 1.08) / dz;
			V = Ez[x_s][y_s][z_v] * dz;
		}
		else if (t > 6.0e-8&&t <= 8.0e-8) {
			Ez[x_s][y_s][z_v] = (t*(1.48 / 2e-8) - 3.6) / dz;
			V = Ez[x_s][y_s][z_v] * dz;
		}
		else if (t > 8.0e-8&&t <= 1.2e-7) {
			Ez[x_s][y_s][z_v] = (t*(0.36 / 4e-8) + 1.6) / dz;
			V = Ez[x_s][y_s][z_v] * dz;
		}
		else if (t > 1.2e-7&&t <= 1.6e-7) {
			Ez[x_s][y_s][z_v] = (t*(0.08 / 4e-8) + 2.44) / dz;
			V = Ez[x_s][y_s][z_v] * dz;
		}
		else {
			Ez[x_s][y_s][z_v] = 2.8 / dz;
			V = 2.8;
		}

		//磁束観測点での電圧
		V1 = 0.0;
		for (k = z_s; k < z_e; k++) {
			V1 += Ez[x_s][y_5][k] * dz;
		}

		//境界条件
		//////////////////////////////////////////// 左面 ////////////////////////////////////////////
		for (i = 0; i <= Nx; i++) {
			for (k = 0; k <= Nz; k++) {
				Ex[i][0][k] = Ex_old[i][1][k] + ((c*dt - dx) / (c*dt + dx))*(Ex[i][1][k] - Ex_old[i][0][k]);
			}
		}
		for (i = 0; i <= Nx; i++) {
			for (k = 0; k <= Nz; k++) {
				Ez[i][0][k] = Ez_old[i][1][k] + ((c*dt - dz) / (c*dt + dz))*(Ez[i][1][k] - Ez_old[i][0][k]);
			}
		}
		//////////////////////////////////////////// 上面 ////////////////////////////////////////////
		for (i = 0; i <= Nx; i++) {
			for (j = 0; j <= Ny; j++) {
				Ex[i][j][Nz] = Ex_old[i][j][Nz - 1] + ((c*dt - dx) / (c*dt + dx))*(Ex[i][j][Nz - 1] - Ex_old[i][j][Nz]);
			}
		}
		for (i = 0; i <= Nx; i++) {
			for (j = 0; j <= Ny; j++) {
				Ey[i][j][Nz] = Ey_old[i][j][Nz - 1] + ((c*dt - dy) / (c*dt + dy))*(Ey[i][j][Nz - 1] - Ey_old[i][j][Nz]);
			}
		}
		//////////////////////////////////////////// 下面 ////////////////////////////////////////////
		for (i = 0; i <= Nx; i++) {
			for (j = 0; j <= Ny; j++) {
				Ex[i][i][0] = Ex_old[i][j][1] + ((c*dt - dy) / (c*dt + dy))*(Ex[i][j][1] - Ex_old[i][j][0]);
			}
		}
		for (i = 0; i <= Nx; i++) {
			for (j = 0; j <= Ny; j++) {
				Ey[i][j][0] = Ey_old[i][j][1] + ((c*dt - dy) / (c*dt + dy))*(Ey[i][j][1] - Ey_old[i][j][0]);
			}
		}
		//////////////////////////////////////////// 手前面 ////////////////////////////////////////////
		for (j = 0; j <= Ny; j++) {
			for (k = 0; k <= Nz; k++) {
				Ey[Nx][j][k] = Ey_old[Nx - 1][j][k] + ((c*dt - dy) / (c*dt + dy))*(Ey[Nx - 1][j][k] - Ey_old[Nx][j][k]);

			}
		}
		for (j = 0; j <= Ny; j++) {
			for (k = 0; k <= Nz; k++) {

				Ez[Nx][j][k] = Ez_old[Nx - 1][j][k] + ((c*dt - dx) / (c*dt + dx))*(Ez[Nx - 1][j][k] - Ez_old[Nx][j][k]);
			}
		}
		//////////////////////////////////////////// 奥面 ////////////////////////////////////////////
		for (j = 0; j <= Ny; j++) {
			for (k = 0; k <= Nz; k++) {
				Ey[0][j][k] = Ey_old[1][j][k] + ((c*dt - dy) / (c*dt + dy))*(Ey[1][j][k] - Ey_old[0][j][k]);

			}
		}
		for (j = 0; j <= Ny; j++) {
			for (k = 0; k <= Nz; k++) {
				Ez[0][j][k] = Ez_old[1][j][k] + ((c*dt - dz) / (c*dt + dz))*(Ez[1][j][k] - Ez_old[0][j][k]);
			}
		}
		//////////////////////////////////////////// 右面 ////////////////////////////////////////////
		for (i = 0; i <= Nx; i++) {
			for (k = 0; k <= Nz; k++) {
				Ex[i][Ny][k] = Ex_old[i][Ny - 1][k] + ((c*dt - dx) / (c*dt + dx))*(Ex[i][Ny - 1][k] - Ex_old[i][Ny][k]);
			}
		}
		for (i = 0; i <= Nx; i++) {
			for (k = 0; k <= Nz; k++) {
				Ez[i][Ny][k] = Ez_old[i][Ny - 1][k] + ((c*dt - dz) / (c*dt + dz))*(Ez[i][Ny - 1][k] - Ez_old[i][Ny][k]);
			}
		}

		t += 0.5*dt;		//時間をdt/2進める

							//磁界の計算//
		for (i = 1; i < Nx; i++) {
			for (j = 0; j < Ny; j++) {
				for (k = 0; k < Nz; k++) {
					Hx[i][j][k] = Hx_old[i][j][k] - (dt / (μ*dy)) * (Ez[i][j + 1][k] - Ez[i][j][k]) + (dt / (μ*dz)) * (Ey[i][j][k + 1] - Ey[i][j][k]);
				}
			}
		}
		for (i = 0; i < Nx; i++) {
			for (j = 1; j < Ny; j++) {
				for (k = 0; k < Nz; k++) {
					Hy[i][j][k] = Hy_old[i][j][k] - (dt / (μ*dz)) * (Ex[i][j][k + 1] - Ex[i][j][k]) + (dt / (μ*dx)) * (Ez[i + 1][j][k] - Ez[i][j][k]);
				}
			}
		}
		for (i = 0; i < Nx; i++) {
			for (j = 0; j < Ny; j++) {
				for (k = 1; k < Nz; k++) {
					Hz[i][j][k] = Hz_old[i][j][k] - (dt / (μ*dx)) * (Ey[i + 1][j][k] - Ey[i][j][k]) + (dt / (μ*dy)) * (Ex[i][j + 1][k] - Ex[i][j][k]);
				}
			}
		}

		//等価媒質定数法(磁界)//
		for (j = y_s + 1; j < y_e; j++) {
			Hx[x_s][j][z_e - 1] = Hx_old[x_s][j][z_e - 1] - (dt / (μm*dy)) * (Ez[x_s][j + 1][z_e - 1] - Ez[x_s][j][z_e - 1]) + (dt / (μm*dz)) * (Ey[x_s][j][z_e] - Ey[x_s][j][z_e - 1]);
			Hx[x_s][j][z_e] = Hx_old[x_s][j][z_e] - (dt / (μm*dy)) * (Ez[x_s][j + 1][z_e] - Ez[x_s][j][z_e]) + (dt / (μm*dz)) * (Ey[x_s][j][z_e + 1] - Ey[x_s][j][z_e]);

			Hz[x_s - 1][j][z_e] = Hz_old[x_s - 1][j][z_e] - (dt / (μm*dx)) * (Ey[x_s][j][z_e] - Ey[x_s - 1][j][z_e]) + (dt / (μm*dy)) * (Ex[x_s - 1][j + 1][z_e] - Ex[x_s - 1][j][z_e]);
			Hz[x_s][j][z_e] = Hz_old[x_s][j][z_e] - (dt / (μm*dx)) * (Ey[x_s + 1][j][z_e] - Ey[x_s][j][z_e]) + (dt / (μm*dy)) * (Ex[x_s][j + 1][z_e] - Ex[x_s][j][z_e]);
		}

		t += 0.5*dt;		//時間をdt/2進める

							//電界・磁界の更新
		for (i = 0; i <= Nx; i++) {
			for (j = 0; j <= Ny; j++) {
				for (k = 0; k <= Nz; k++) {
					Ex_old[i][j][k] = Ex[i][j][k];
					Ey_old[i][j][k] = Ey[i][j][k];
					Ez_old[i][j][k] = Ez[i][j][k];
					Hx_old[i][j][k] = Hx[i][j][k];
					Hy_old[i][j][k] = Hy[i][j][k];
					Hz_old[i][j][k] = Hz[i][j][k];
				}
			}
		}

		//Vf = Ex[x_s][y_e][z_lr] * ds;		//電圧
		I = (Hx[x_s][y_e - 1][z_lr] - Hx[x_s][y_e][z_lr])*dx + (Hy[x_s][y_e][z_lr] - Hy[x_s - 1][y_e][z_lr])*dy;		//電流
		H = Hx[x_s][y_5][z_v + 1];


		printf("時間");
		printf("%e\t\t\t", t);
		printf("電流：");
		printf("%e\t\t", I);
		printf("電圧：");
		printf("%e\t\t", V);
		printf("インピーダンス：%f\n\n", V / I);

		fprintf(fp, "%e,%e\n", t, V);
		fprintf(fp1, "%e,%e\n", t, I);
		fprintf(fp2, "%e,%e\n", t, V1);
		fprintf(fp3, "%e, %e\n", t, V1);
		fprintf(fp4, "%e %e\n", t, I);
		fprintf(fp5, "%e %e\n", t, H);


	}

	fclose(fp);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);

	T2 = clock();
	T3 = T2 - T1;
	printf("%d\n\n\n", T3 / CLOCKS_PER_SEC);

}