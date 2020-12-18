#include<stdio.h>
#include<math.h>
#include<time.h>

//���[�v��(x,y,z)//
int Nx = 40;
int Ny = 251;
int Nz = 40;

/*��͗̈�*/
//x=120[mm],y=753[mm],z=120[mm]

//x�����̓���
#define x_s  20			//x��������
//y�����̓���
#define y_s  17			//y���������n�܂�(���[����51[mm]�̏ꏊ)
#define y_5  164		//y����(���E����ꏊ�@���[����51�{441[mm])
#define y_e  234		//y���������I���(�E�[����124[mm]�̏ꏊ)
//z�����̓���
#define z_s  17			//z���������n�܂�
#define z_v  20			//z���������`�d��
#define z_e  23			//z���������I���
#define z_r  17			//z���������`��R
#define z_re 23			//z������R�`����
//��R�̈ʒu
#define z_lr 21			//��R�̏ꏊ
//�C���_�N�^���X�̈ʒu
#define z_ll 19			//�C���_�N�^���X�̏ꏊ

//��//
double pi = 4.0 * atan(1.0);

//�T�C�Y//
double dx = 3.0e-3;
double dy = 3.0e-3;
double dz = 3.0e-3;
double dt = 5.0e-12;
double ds = 3.0e-3;

//�p�����[�^�̒l//
double c = 3.0e8;					//�d���g�̐i�s���x�@c[m/sec]
double �� = 1.0e-8;					//���d���@��[S/m]
double �� = 4.0e-7*pi;				//�������@��[H/m]
double �� = 1.0e-9 / (36 * pi);		//�U�d���@��[F/m]
double t;
double R = 1370;					//��R	R[��]

double L = 5.94e-3;					//�C���_�N�^���X L[H]

									//�d�E,���E//
double Ex[41][252][41], Ex_old[41][252][41];
double Ey[41][252][41], Ey_old[41][252][41];
double Ez[41][252][41], Ez_old[41][252][41];

double Hx[41][252][41], Hx_old[41][252][41];
double Hy[41][252][41], Hy_old[41][252][41];
double Hz[41][252][41], Hz_old[41][252][41];

//�d��,�d��
double V = 0.0;
double V1 = 0.0;
double Vf = 0.0;
double I = 0.0;
double H = 0.0;

//�����̐ݒ�
double r = 1.0e-3;			//�����̔��a[m]
double r0 = 0.230;			//�����̓������a[m]
double d = 1.8e-2;			//���s2���̋���[m]

							//�����}���萔�@
double m = 0;				//�����}���萔
double ��m = 0;				//�����ߖT�ł̃�m[F/m]
double ��m = 0;				//�����ߖT�ł̃�m[H/m]


							//�d�E�̌v�Z�̌W��//
double E1 = 0;		double E2 = 0;
//�����}���萔�@��K���������̓d�E�̌v�Z�̌W��//
double Em1 = 0;		double Em2 = 0;

//��R�̓d�E�̌v�Z�̌W��//
double R1 = 0;      double R2 = 0;

//�C���_�N�^���X�̓d�E�̌v�Z�̌W��//
double integral = 0;
//�C���_�N�^���X�̗��_�l//
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

	sprintf_s(moji, "�d���̎��ԕω�(im).csv", count);
	fopen_s(&fp, moji, "w");

	sprintf_s(moji, "�d���̎��ԕω�(im).csv", count);
	fopen_s(&fp1, moji, "w");

	sprintf_s(moji, "���E�ϑ��_�̏ꏊ�̓d��(im).csv", count);
	fopen_s(&fp2, moji, "w");

	sprintf_s(moji, "���E�ϑ��_�̏ꏊ�̓d��(im).txt", count);
	fopen_s(&fp3, moji, "w");

	sprintf_s(moji, "�d���̎��ԕω�(im).txt", count);
	fopen_s(&fp4, moji, "w");

	sprintf_s(moji, "���E�̎��ԕω�(im).txt", count);
	fopen_s(&fp5, moji, "w");

	m = (log(1 / r0) / log(ds / r));		//�����}���萔
	��m = ��*m;							    //�����ߖT�ł̃�m[F/m]
	��m = �� / m;						    //�����ߖT�ł̃�m[H/m]

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

	E1 = (2 * �� - ��*dt) / (2 * �� + ��*dt);			E2 = (2 * dt) / (��*dt + 2 * ��);			//�d�E�̌v�Z�̌W��//
	Em1 = (2 * ��m - ��*dt) / (2 * ��m + ��*dt);		Em2 = (2 * dt) / (��*dt + 2 * ��m);			//�����}���萔�@��K���������̓d�E�̌v�Z�̌W��//


	R1 = (2 * R*��*dx*dy - dt * dz) / (2 * R*��*dx*dy + dt * dz);
	R2 = dt*(2 * R*dx*dy) / (2 * ��*R*dx*dy + dz*dt);						//��R�̓d�E�̌v�Z�̌W��//

	for (t = 0; t <= 1.0e-6; ) {

		//�d�E�̌v�Z//
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

		//��R�̓d�E�v�Z//
		Ez[x_s][y_e][z_lr] = R1*Ez_old[x_s][y_e][z_lr] + (R2 / dx)*(Hy_old[x_s][y_e][z_lr] - Hy_old[x_s - 1][y_e][z_lr]) - (R2 / dy)*(Hx_old[x_s][y_e][z_lr] - Hx_old[x_s][y_e - 1][z_lr]);

		//�C���_�N�^���X�̓d�E�v�Z//
		Ez[x_s][y_e][z_ll] = 1.0*Ez_old[x_s][y_e][z_ll] + (dt / ��)*1.0*(Hy_old[x_s][y_e][z_ll] - Hy_old[x_s - 1][y_e][z_ll]) - (dt / ��)*1.0*(Hx_old[x_s][y_e][z_ll] - Hx_old[x_s][y_e - 1][z_ll]) - integral;

		//�C���_�N�^���X�̌v�Z�̌W��//
		integral += dz*dt*dt / (��*L*dx*dy)*Ez[x_s][y_e][z_ll];

		//�����}���萔�@(�d�E)//
		for (j = y_s + 1; j < y_e; j++) {
			Ex[x_s - 1][j][z_e] = Em1 * Ex_old[x_s - 1][j][z_e] + (Em2 / dy) * (Hz_old[x_s - 1][j][z_e] - Hz_old[x_s - 1][j - 1][z_e]) - (Em2 / dz) * (Hy_old[x_s - 1][j][z_e] - Hy_old[x_s - 1][j][z_e - 1]);
			Ex[x_s][j][z_e] = Em1 * Ex_old[x_s][j][z_e] + (Em2 / dy) * (Hz_old[x_s][j][z_e] - Hz_old[x_s][j - 1][z_e]) - (Em2 / dz) * (Hy_old[x_s][j][z_e] - Hy_old[x_s][j][z_e - 1]);

			Ez[x_s][j][z_e - 1] = Em1 * Ez_old[x_s][j][z_e - 1] + (Em2 / dx) * (Hy_old[x_s][j][z_e - 1] - Hy_old[x_s - 1][j][z_e - 1]) - (Em2 / dy) * (Hx_old[x_s][j][z_e - 1] - Hx_old[x_s][j - 1][z_e - 1]);
			Ez[x_s][j][z_e] = Em1 * Ez_old[x_s][j][z_e] + (Em2 / dx) * (Hy_old[x_s][j][z_e] - Hy_old[x_s - 1][j][z_e]) - (Em2 / dy) * (Hx_old[x_s][j][z_e] - Hx_old[x_s][j - 1][z_e]);
		}

		//////////////////////////////////////////// ����������(��̓���) ////////////////////////////////////////////

		for (j = y_s; j < y_e; j++) {
			Ey[x_s][j][z_e] = 0.0;
		}

		//////////////////////////////////////////// ����������(���̓���) ////////////////////////////////////////////

		for (j = y_s; j < y_e; j++) {
			Ey[x_s][j][z_s] = 0.0;
		}

		//////////////////////////////////////////// ����z���� ////////////////////////////////////////////

		for (k = z_s; k < z_v; k++) {			//���̓����`�d��
			Ez[x_s][y_s][k] = 0.0;
		}
		for (k = z_v + 1; k < z_e; k++) { 		//�d���`��̓���
			Ez[x_s][y_s][k] = 0.0;
		}
		for (k = z_r; k < z_ll; k++) {			//���̓����`�R�C��
			Ez[x_s][y_e][k] = 0.0;
		}
		for (k = z_ll + 1; k < z_lr; k++) {		//�R�C���`��R
			Ez[x_s][y_e][k] = 0.0;
		}
		for (k = z_lr + 1; k < z_re; k++) {      //�R�C���`����
			Ez[x_s][y_e][k] = 0.0;
		}

		//�d�E����d����^���鎮
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

		//�����ϑ��_�ł̓d��
		V1 = 0.0;
		for (k = z_s; k < z_e; k++) {
			V1 += Ez[x_s][y_5][k] * dz;
		}

		//���E����
		//////////////////////////////////////////// ���� ////////////////////////////////////////////
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
		//////////////////////////////////////////// ��� ////////////////////////////////////////////
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
		//////////////////////////////////////////// ���� ////////////////////////////////////////////
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
		//////////////////////////////////////////// ��O�� ////////////////////////////////////////////
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
		//////////////////////////////////////////// ���� ////////////////////////////////////////////
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
		//////////////////////////////////////////// �E�� ////////////////////////////////////////////
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

		t += 0.5*dt;		//���Ԃ�dt/2�i�߂�

							//���E�̌v�Z//
		for (i = 1; i < Nx; i++) {
			for (j = 0; j < Ny; j++) {
				for (k = 0; k < Nz; k++) {
					Hx[i][j][k] = Hx_old[i][j][k] - (dt / (��*dy)) * (Ez[i][j + 1][k] - Ez[i][j][k]) + (dt / (��*dz)) * (Ey[i][j][k + 1] - Ey[i][j][k]);
				}
			}
		}
		for (i = 0; i < Nx; i++) {
			for (j = 1; j < Ny; j++) {
				for (k = 0; k < Nz; k++) {
					Hy[i][j][k] = Hy_old[i][j][k] - (dt / (��*dz)) * (Ex[i][j][k + 1] - Ex[i][j][k]) + (dt / (��*dx)) * (Ez[i + 1][j][k] - Ez[i][j][k]);
				}
			}
		}
		for (i = 0; i < Nx; i++) {
			for (j = 0; j < Ny; j++) {
				for (k = 1; k < Nz; k++) {
					Hz[i][j][k] = Hz_old[i][j][k] - (dt / (��*dx)) * (Ey[i + 1][j][k] - Ey[i][j][k]) + (dt / (��*dy)) * (Ex[i][j + 1][k] - Ex[i][j][k]);
				}
			}
		}

		//�����}���萔�@(���E)//
		for (j = y_s + 1; j < y_e; j++) {
			Hx[x_s][j][z_e - 1] = Hx_old[x_s][j][z_e - 1] - (dt / (��m*dy)) * (Ez[x_s][j + 1][z_e - 1] - Ez[x_s][j][z_e - 1]) + (dt / (��m*dz)) * (Ey[x_s][j][z_e] - Ey[x_s][j][z_e - 1]);
			Hx[x_s][j][z_e] = Hx_old[x_s][j][z_e] - (dt / (��m*dy)) * (Ez[x_s][j + 1][z_e] - Ez[x_s][j][z_e]) + (dt / (��m*dz)) * (Ey[x_s][j][z_e + 1] - Ey[x_s][j][z_e]);

			Hz[x_s - 1][j][z_e] = Hz_old[x_s - 1][j][z_e] - (dt / (��m*dx)) * (Ey[x_s][j][z_e] - Ey[x_s - 1][j][z_e]) + (dt / (��m*dy)) * (Ex[x_s - 1][j + 1][z_e] - Ex[x_s - 1][j][z_e]);
			Hz[x_s][j][z_e] = Hz_old[x_s][j][z_e] - (dt / (��m*dx)) * (Ey[x_s + 1][j][z_e] - Ey[x_s][j][z_e]) + (dt / (��m*dy)) * (Ex[x_s][j + 1][z_e] - Ex[x_s][j][z_e]);
		}

		t += 0.5*dt;		//���Ԃ�dt/2�i�߂�

							//�d�E�E���E�̍X�V
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

		//Vf = Ex[x_s][y_e][z_lr] * ds;		//�d��
		I = (Hx[x_s][y_e - 1][z_lr] - Hx[x_s][y_e][z_lr])*dx + (Hy[x_s][y_e][z_lr] - Hy[x_s - 1][y_e][z_lr])*dy;		//�d��
		H = Hx[x_s][y_5][z_v + 1];


		printf("����");
		printf("%e\t\t\t", t);
		printf("�d���F");
		printf("%e\t\t", I);
		printf("�d���F");
		printf("%e\t\t", V);
		printf("�C���s�[�_���X�F%f\n\n", V / I);

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