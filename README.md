# Bridge-Stress-Analysis-and-Evaluation-Program
# This program is Python-based and calculates bridge stress distribution, incorporates crack data, and evaluates failure using Mohr’s Criterion. Key features include stress tensor computation, an intuitive GUI, and visualizations such as stress diagrams, crack locations, and safety zones, effectively supporting bridge safety assessment.
import matplotlib
matplotlib.use('TkAgg')  # GUI 백엔드 설정

import matplotlib.pyplot as plt
import sympy as sp
import numpy as np

# 한글 폰트 설정 및 Unicode Minus 설정
plt.rcParams['font.family'] = 'Malgun Gothic'  # 윈도우의 경우 'Malgun Gothic' 사용
plt.rcParams['axes.unicode_minus'] = False     # 유니코드 마이너스 대신 일반 하이픈 사용

# 교량 길이 입력
bridge_length = float(input("교량 길이(m)를 입력하세요: "))

# 단면의 폭과 높이 입력
width = float(input("단면의 폭(m)를 입력하세요: "))
height = float(input("단면의 높이(m)를 입력하세요: "))

# 하중 정보 입력
loads = []
while True:
    load_type = input("하중 타입을 입력하세요 (point_moment/point_load/distributed_load, 입력기준(moment : 반시계 기준 입력, load : 아래방향 표시된 값 입력) 종료하려면 '-1' 입력): ")
    if load_type == '-1':
        break
    elif load_type == 'point_moment':
        pos = float(input("포인트 모멘트 하중의 x축 위치(m)를 입력하세요: "))
        magnitude = float(input("하중의 크기(kN*m): "))
        loads.append({'type': 'point_moment', 'position': pos, 'magnitude': magnitude})
    elif load_type == 'point_load':
        pos = float(input("포인트 수직 하중의 x축 위치(m)를 입력하세요: "))
        magnitude = float(input("하중의 크기(kN): "))
        loads.append({'type': 'point_load', 'position': pos, 'magnitude': magnitude})
    elif load_type == 'distributed_load':
        start_pos = float(input("연속 수직 하중의 시작 x축 위치(m)를 입력하세요: "))
        end_pos = float(input("연속 수직 하중의 끝 x축 위치(m)를 입력하세요: "))
        magnitude = float(input("하중의 크기(kN/m): "))
        loads.append({'type': 'distributed_load', 'start': start_pos, 'end': end_pos, 'magnitude': magnitude})

# SymPy를 사용하여 반력 계산
R_A, R_B = sp.symbols('R_A R_B')
moment_eq = 0  # 모멘트 방정식 초기화
force_eq = R_A + R_B  # 수직력 평형 초기화

# 하중에 따른 수직력과 모멘트 방정식 구성
for load in loads:
    if load['type'] == 'point_load':
        force_eq -= load['magnitude']
        moment_eq -= load['magnitude'] * (bridge_length - load['position'])
    elif load['type'] == 'distributed_load':
        dist_force = load['magnitude'] * (load['end'] - load['start'])
        centroid = (load['start'] + load['end']) / 2
        force_eq -= dist_force
        moment_eq -= dist_force * (bridge_length - centroid)
    elif load['type'] == 'point_moment':
        moment_eq -= load['magnitude']

# B 지점(bridge_length 위치)에서의 모멘트 평형: 회전평형 사용
moment_eq = sp.Eq(moment_eq + R_A * bridge_length, 0)
force_eq = sp.Eq(force_eq, 0)

# 연립 방정식 풀기
solution = sp.solve((moment_eq, force_eq), (R_A, R_B))
R_A_value = solution[R_A]
R_B_value = solution[R_B]

print(f"반력 R_A: {R_A_value:.2f} kN")
print(f"반력 R_B: {R_B_value:.2f} kN")

# 자유물체도(FBD) 그리기
fig, ax = plt.subplots(figsize=(10, 4))

# 빔 그리기
ax.plot([0, bridge_length], [0, 0], 'k', linewidth=5)

# 반력 표시 (위 방향 화살표로 RA, RB)
R_A_pos = 0
R_B_pos = bridge_length

ax.arrow(R_A_pos, -1, 0, 0.8, head_width=0.2, head_length=0.5, fc='b', ec='b')
ax.text(R_A_pos, -1.5, f'R_A = {R_A_value:.2f} kN', color='b', ha='center')

ax.arrow(R_B_pos, -1, 0, 0.8, head_width=0.2, head_length=0.5, fc='r', ec='r')
ax.text(R_B_pos, -1.5, f'R_B = {R_B_value:.2f} kN', color='r', ha='center')

# 하중 표시
for load in loads:
    if load['type'] == 'point_load':
        ax.arrow(load['position'], -1, 0, 0.8, head_width=0.2, head_length=0.5, fc='g', ec='g')
        ax.text(load['position'], -1.5, f'-{load["magnitude"]} kN', color='g', ha='center')
    elif load['type'] == 'distributed_load':
        num_arrows = 5
        step = (load['end'] - load['start']) / (num_arrows - 1)
        for i in range(num_arrows):
            pos = load['start'] + i * step
            ax.arrow(pos, -0.3, 0, 0.2, head_width=0.1, head_length=0.2, fc='g', ec='g')
        ax.plot([load['start'], load['end']], [-0.3, -0.3], 'g--', linewidth=1)
        ax.text((load['start'] + load['end']) / 2, -0.5, f'-{load["magnitude"]} kN/m', color='g', ha='center')
    elif load['type'] == 'point_moment':
        circle = plt.Circle((load['position'], 0), 0.5, color='purple', fill=False, linestyle='--')
        ax.add_artist(circle)
        ax.plot(load['position'], 0.5, color='purple', marker=(3, 0, -90), markersize=15)
        ax.text(load['position'], 0.7, f'-{load["magnitude"]} kN*m', color='purple', ha='center')

ax.set_ylim(-2, 2)
ax.set_xlim(0, bridge_length)
ax.axis('off')
plt.title("Free-Body Diagram with External Sign Convention (kN, m)")
plt.show()

import sympy as sp

# 변수 정의 및 초기화
x = sp.symbols('x')
w_x = 0  # 초기 하중 함수 설정

# 사용자 입력에 따른 하중 정의 - DiracDelta와 Heaviside를 사용
moments = []  # 별도의 point_moment 저장용 리스트 추가
for load in loads:
    if load['type'] == 'point_load':
        # 포인트 하중: DiracDelta로 표현
        w_x -= load['magnitude'] * sp.DiracDelta(x - load['position'])
    elif load['type'] == 'distributed_load':
        # 분포 하중: Heaviside로 구간 설정
        w_x -= load['magnitude'] * (sp.Heaviside(x - load['start']) - sp.Heaviside(x - load['end']))
    elif load['type'] == 'point_moment':
        # 모멘트는 전단력/모멘트 방정식에서 직접 처리하므로 따로 저장
        moments.append({'position': load['position'], 'magnitude': load['magnitude']})

# 반력 점 하중 추가: DiracDelta 함수로 특정 위치에서만 작용하도록 설정
w_x += R_A_value * sp.DiracDelta(x) + R_B_value * sp.DiracDelta(x - bridge_length)

# 전단력 V(x) 계산
V_x = sp.integrate(w_x, x) 

# 모멘트 M(x) 계산
M_x = sp.integrate(V_x, x) 

# point_moment 추가 반영: 모멘트 함수에서 순간적인 불연속 반영
for moment in moments:
    M_x += moment['magnitude'] * sp.Heaviside(x - moment['position'])


# 각 함수 간소화
w_x = sp.simplify(w_x)
V_x = sp.simplify(V_x)
M_x = sp.simplify(M_x)

# 함수 확인용 출력
print("전체 하중 함수 w(x):")
sp.pprint(w_x)
print("\n전체 전단력 함수 V(x):")
sp.pprint(V_x)
# print("\n전체 모멘트 함수 M(x):")
# sp.pprint(M_x) 깔끔하게 만들기 위해 안보이도록 만듬

# 기존 모멘트와 전단력 계산
V_x = sp.integrate(w_x, x)
M_x = sp.integrate(V_x, x)
V_x = sp.simplify(V_x)
M_x = sp.simplify(M_x)

# 기존의 콘솔 출력 부분 제거
# print("전체 모멘트 함수 M(x):")
# sp.pprint(M_x)

# GUI 코드 추가
import tkinter as tk
from tkinter import messagebox

def show_moment_equation():
    """
    버튼 클릭 시 모멘트 식을 팝업 창으로 보여줍니다.
    """
    global M_x
    if M_x is not None:
        moment_equation = sp.pretty(M_x)  # SymPy 식을 읽기 쉽게 변환
        messagebox.showinfo("모멘트 식", f"모멘트 M(x):\n\n{moment_equation}")
    else:
        messagebox.showwarning("경고", "모멘트 식이 아직 계산되지 않았습니다!")

# Tkinter GUI 초기화
root = tk.Tk()
root.title("교량 상태 평가 프로그램")

# 메인 프레임 생성
frame = tk.Frame(root)
frame.pack(pady=10, padx=10)

# 제목 라벨
title_label = tk.Label(frame, text="교량 상태 평가 프로그램", font=("Arial", 16, "bold"))
title_label.pack(pady=5)

# 모멘트 식 보기 버튼
moment_button = tk.Button(frame, text="모멘트 식 보기", command=show_moment_equation, font=("Arial", 12))
moment_button.pack(pady=5)

# 종료 버튼
exit_button = tk.Button(frame, text="종료", command=root.quit, font=("Arial", 12))
exit_button.pack(pady=5)

# Tkinter 메인 루프 실행
root.mainloop()

def get_valid_z_position(width):
    """
    사용자가 입력한 z 위치가 단면 폭 범위를 초과할 경우 다시 입력받습니다.

    Parameters:
        width (float): 단면 폭 (m)

    Returns:
        float: 유효한 z 위치 값
    """
    while True:
        try:
            z_position = float(input(f"응력을 계산할 위치 z (m)를 입력하세요 (범위: {-width/2} ~ {width/2}): "))
            if -width / 2 <= z_position <= width / 2:
                return z_position
            else:
                print(f"오류: 입력된 z 위치가 단면 폭 범위를 초과했습니다. 범위는 {-width/2} ~ {width/2}입니다.")
        except ValueError:
            print("잘못된 입력입니다. 숫자를 입력해주세요.")
            
# 임의의 위치 (x, y, z) 입력받기
x_val = float(input("응력을 계산할 위치 x (m)를 입력하세요: "))
y_val = float(input("응력을 계산할 위치 y (m)를 입력하세요: "))
z_val = get_valid_z_position(width)  # 유효한 z 위치를 반환받음

# 단면 2차 모멘트 계산 (I = width * height^3 / 12)
I = width * height**3 / 12
print(f"단면 2차 모멘트 I: {I:.4f} mm^4")

# 응력 텐서 계산 함수 정의
def compute_stress_tensor(moment, shear_force, height, width, y_position):
    """
    주어진 조건에서 3x3 응력 텐서를 계산합니다.

    Parameters:
        moment (float): 모멘트 (kN·m)
        shear_force (float): 전단력 (kN)
        height (float): 단면 높이 (m)
        width (float): 단면 폭 (m)
        y_position (float): y 위치 (m)
    
    Returns:
        numpy.ndarray: 3x3 응력 텐서
    """
    # y 위치가 단면 높이의 절반을 초과하면 오류 처리
    if abs(y_position) > height / 2:
        print("오류: 입력된 위치 y가 단면 범위를 초과했습니다. (y는 단면 높이의 절반 이하여야 합니다)")
        return y_val

    # 단면 2차 모멘트 (I) 계산
    inertia_moment = (width * height**3) / 12

    # 수직 응력 계산
    axial_stress = -(moment * y_position) / inertia_moment

    # 전단 응력을 계산하기 위한 단면 1차 모멘트 (Q)와 두께 설정
    thickness = width
    first_moment_area = (width / 2) * ((height**2 / 4) - y_position**2)
    
    # 전단 응력 계산
    shear_stress = -(shear_force * first_moment_area) / (inertia_moment * thickness)

    # 3x3 응력 텐서 생성
    stress_tensor = np.array([
        [axial_stress, shear_stress, 0],  # x 방향
        [shear_stress, 0, 0],            # y 방향
        [0, 0, 0]                        # z 방향
    ])
    return stress_tensor

# 응력 텐서 계산 및 출력
moment_at_x = M_x.subs(x, x_val)
shear_force_at_x = V_x.subs(x, x_val)
stress_tensor = compute_stress_tensor(moment_at_x, shear_force_at_x, height, width, y_val)

# SymPy 객체를 숫자 배열로 변환하여 소수점 둘째 자리까지 출력
if stress_tensor is not None:
    print("\n위치 (x={}, y={})에서의 응력 텐서(kPa):".format(x_val, y_val))
    for row in stress_tensor:
        formatted_row = [f"{float(value):.2f}" for value in row]
        print("[", " ".join(formatted_row), "]")

# 프로그램 종료 방지
input("Press Enter to exit...")

# lambdify를 사용하여 numpy가 아닌 sympy 모듈로 변환
V_func = sp.lambdify(x, V_x, modules='sympy')
M_func = sp.lambdify(x, M_x, modules='sympy')

# x 값 범위 설정 및 계산
x_vals = np.linspace(0, bridge_length, 1000)
V_vals = np.array([float(V_func(val)) for val in x_vals])
M_vals = np.array([float(M_func(val)) for val in x_vals])

# 전단력 선도
plt.figure(figsize=(10, 4))
plt.plot(x_vals, V_vals, label="전단력 V(x)", color='green')
plt.xlabel("x (m)")
plt.ylabel("전단력 V(x) (kN)")
plt.title("전단력 선도")
plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
plt.xlim(0, bridge_length)
plt.legend()
plt.grid(True)
plt.show()

# 모멘트 선도
plt.figure(figsize=(10, 4))
plt.plot(x_vals, M_vals, label="모멘트 M(x)", color='orange')
plt.xlabel("x (m)")
plt.ylabel("모멘트 M(x) (kN·m)")
plt.title("모멘트 선도")
plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
plt.xlim(0, bridge_length)
plt.legend()
plt.grid(True)
plt.show()

#응력분포도 플롯
# 필요 변수 및 함수 정의
y_vals = np.linspace(-height/2, height/2, 100)  # 단면 높이를 따라 y 위치를 정의

# 각 y 위치에서 응력을 계산하는 함수
def compute_stress_distribution(moment, shear_force, height, width):
    axial_stresses = []
    shear_stresses = []
    
    for y_position in y_vals:
        # 수직 응력 계산
        axial_stress = -(moment * y_position) / I
        
        # 전단 응력 계산
        first_moment_area = (width / 2) * ((height**2 / 4) - y_position**2)
        shear_stress = -(shear_force * first_moment_area) / (I * width)
        
        axial_stresses.append(axial_stress)
        shear_stresses.append(shear_stress)
    
    return np.array(axial_stresses), np.array(shear_stresses)

# x 위치에서의 모멘트와 전단력
moment_at_x = M_x.subs(x, x_val)
shear_force_at_x = V_x.subs(x, x_val)

# x와 y의 격자 생성
x_vals = np.linspace(-width / 2, width / 2, 100)
y_vals = np.linspace(-height / 2, height / 2, 100)
X, Y = np.meshgrid(x_vals, y_vals)

# 응력 계산 함수
def compute_stress_field(moment, shear_force, X, Y, I, width):
    # 축방향 응력 σ_xx
    axial_stress = -(moment * Y) / I
    # 전단 응력 τ_xy
    first_moment_area = (width / 2) * ((height ** 2 / 4) - Y ** 2)
    shear_stress = -(shear_force * first_moment_area) / (I * width)
    return axial_stress, shear_stress

# x 위치에서의 모멘트와 전단력
moment_at_x = float(M_x.subs(x, x_val))
shear_force_at_x = float(V_x.subs(x, x_val))

# 축방향 응력과 전단 응력 계산
axial_stress_field, shear_stress_field = compute_stress_field(
    moment_at_x, shear_force_at_x, X, Y, I, width
)

# 축방향 응력과 전단 응력 계산 및 입력 위치 표시
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# 축방향 응력 분포 (σ_xx)
contour1 = ax[0].contourf(X, Y, axial_stress_field, cmap='coolwarm', levels=50)
ax[0].set_title("축방향 응력 (σ_xx)")
ax[0].set_xlabel("z-단면의 폭(mm)")
ax[0].set_ylabel("y-단면의 높이(mm)")
ax[0].axhline(0, color='gray', linestyle='--', linewidth=0.5)
ax[0].axvline(0, color='gray', linestyle='--', linewidth=0.5)
# 입력 위치 표시
ax[0].scatter(z_val, y_val, color='black', s=50, label=f"입력 위치: z={z_val:.2f}, y={y_val:.2f}")
ax[0].legend()
fig.colorbar(contour1, ax=ax[0], label="응력 (kPa)")

# 전단 응력 분포 (τ_xy)
contour2 = ax[1].contourf(X, Y, shear_stress_field, cmap='coolwarm', levels=50)
ax[1].set_title("전단 응력 (τ_xy)")
ax[1].set_xlabel("z-단면의 폭(mm)")
ax[1].set_ylabel("y-단면의 높이(mm)")
ax[1].axhline(0, color='gray', linestyle='--', linewidth=0.5)
ax[1].axvline(0, color='gray', linestyle='--', linewidth=0.5)
# 입력 위치 표시
ax[1].scatter(z_val, y_val, color='black', s=50, label=f"입력 위치: z={z_val:.2f}, y={y_val:.2f}")
ax[1].legend()
fig.colorbar(contour2, ax=ax[1], label="응력 (kPa)")

plt.tight_layout()
plt.show()

#응력텐서를 eigenvalue를 이용하여 주응력 구하기
import numpy as np

# 주어진 응력 텐서에서 주응력을 계산하는 함수
import numpy as np

# 주어진 응력 텐서에서 주응력을 계산하는 함수
# 주응력 계산 함수 수정
def calculate_principal_stresses_2d_from_3x3(stress_tensor):
    """
    3x3 응력 텐서에서 2x2 성분만 사용하여 주응력을 계산합니다.
    z-방향 응력(σ3)은 항상 0으로 가정합니다.

    Parameters:
        stress_tensor (numpy.ndarray): 3x3 응력 텐서

    Returns:
        tuple: (sigma1, sigma2, sigma3) 주응력
    """
    # 2x2 응력 텐서 추출
    stress_tensor_2d = stress_tensor[:2, :2]

    # SymPy 객체를 부동소수점(float)으로 변환
    stress_tensor_2d = np.array(stress_tensor_2d, dtype=float)

    # 2x2 응력 텐서의 고유값 계산
    eigenvalues = np.linalg.eigvalsh(stress_tensor_2d)

    # 고유값 내림차순 정렬
    sigma1, sigma2 = sorted(eigenvalues, reverse=True)

    # z축 주응력은 0으로 설정
    sigma3 = 0

    return sigma1, sigma2, sigma3

# 응력 텐서 계산 및 출력
moment_at_x = M_x.subs(x, x_val)
shear_force_at_x = V_x.subs(x, x_val)
stress_tensor = compute_stress_tensor(moment_at_x, shear_force_at_x, height, width, y_val)

# SymPy 객체를 숫자 배열로 변환하여 소수점 둘째 자리까지 출력
if stress_tensor is not None:
    # SymPy 객체를 float로 변환
    stress_tensor = np.array(stress_tensor, dtype=float)

    print("\n위치 (x={}, y={})에서의 응력 텐서(kPa):".format(x_val, y_val))
    for row in stress_tensor:
        formatted_row = [f"{value:.2f}" for value in row]
        print("[", " ".join(formatted_row), "]")
    
    # 주응력 계산 (2D 최적화)
    sigma1, sigma2, sigma3 = calculate_principal_stresses_2d_from_3x3(stress_tensor)

    # 검증: σ1 >= σ2 확인
    if sigma1 < sigma2:
        raise ValueError(f"오류: 계산된 주응력이 σ1 >= σ2 조건을 만족하지 않습니다! (σ1={sigma1}, σ2={sigma2})")

    # 결과 출력
    #(가장 큰 축 방향 응력)
    print(f"\n주응력 σ1(kPa): {sigma1:.2f} ")
    #(다음으로 큰 응력)
    print(f"주응력 σ2(kPa): {sigma2:.2f} ")
    #(z-방향 응력)
    print(f"주응력 σ3(kPa): {sigma3:.2f} ")

# 크랙 정보를 입력받기
cracks = []
while True:
    crack_input = input("크랙 정보를 입력하세요 (x 위치, 하단에서 깊이 d, 종료하려면 -1 입력): ")
    if crack_input == '-1':
        break
    try:
        x_position, depth = map(float, crack_input.split(','))
        cracks.append({'x': x_position, 'd': depth})
    except ValueError:
        print("잘못된 입력입니다. x 위치와 깊이를 쉼표로 구분해 입력하세요.")

# 크랙의 영향을 고려한 단면 2차 모멘트 및 중립축 계산 함수
def compute_effective_properties_with_crack(width, height, crack_depth):
    """
    크랙을 고려한 유효 단면 높이, 중립축, 2차 모멘트 계산
    
    Parameters:
        width (float): 단면 폭
        height (float): 단면 높이
        crack_depth (float): 크랙 깊이

    Returns:
        h_eff (float): 유효 단면 높이
        neutral_axis (float): 중립축 위치 (글로벌 좌표계)
        I_eff (float): 유효 단면 2차 모멘트
    """
    h_eff = height - crack_depth  # 남은 유효 높이
    A_eff = width * h_eff  # 유효 단면 면적
    y_c = -(height / 2) + (h_eff / 2)  # 글로벌 좌표계에서 중립축 위치
    
    # 유효 2차 모멘트 계산 (중립축 이동 고려)
    I_eff = (width * h_eff**3) / 12 + A_eff * (abs(crack_depth - y_c)**2)
    
    return h_eff, y_c, I_eff

# 유효 단면 2차 모멘트 출력 및 크랙을 고려한 응력 텐서 계산 함수
def compute_stress_tensor_with_cracks(moment, shear_force, width, height, y_position, cracks, x_position):
    """
    크랙 정보를 고려하여 응력 텐서를 계산합니다. 크랙이 없는 경우 기존 단면 기준으로 계산합니다.
    """
    # x 위치에서 가장 가까운 크랙 확인
    crack = next((c for c in cracks if abs(c['x'] - x_position) < 1e-6), None)
    
    if crack is None:
        # 크랙이 없는 경우 기본 단면 기준으로 기존 응력 텐서 계산
        print(f"x = {x_position} 위치에서 크랙이 없습니다. 기본 단면 기준으로 응력을 계산합니다.")
        return compute_stress_tensor(moment, shear_force, height, width, y_position)
    else:
        # 크랙이 있는 경우 유효 단면 계산
        crack_depth = crack['d']
        h_eff, neutral_axis, I_eff = compute_effective_properties_with_crack(width, height, crack_depth)
        
        # 크랙 내부인지 확인
        crack_limit = -height / 2 + crack_depth  # 크랙 깊이에 대한 기준
        if y_position <= crack_limit:
            print(f"y = {y_position}가 크랙 내부에 위치합니다. 응력은 0으로 계산됩니다.")
            return np.zeros((3, 3))
    
    # 유효 중립축 기준으로 y_local 계산
    y_local = y_position - neutral_axis
    
    # 축방향 응력 계산
    axial_stress = -(moment * y_local) / I_eff

    # 전단 응력을 계산하기 위한 단면 1차 모멘트 (Q)
    Q = (width / 2) * ((h_eff**2 / 4) - y_local**2)  # <-- 여기에서 h_eff를 사용
    shear_stress = -(shear_force * Q) / (I_eff * width)

    # 3x3 응력 텐서 작성
    stress_tensor = np.array([
        [axial_stress, shear_stress, 0],
        [shear_stress, 0, 0],
        [0, 0, 0]
    ])
    return stress_tensor

# 응력 텐서 계산 예시
x_val = float(input("응력을 계산할 위치 x (m)를 입력하세요: "))
y_val = float(input("응력을 계산할 위치 y (m)를 입력하세요: "))

moment_at_x = M_x.subs(x, x_val)
shear_force_at_x = V_x.subs(x, x_val)

# 크랙을 고려한 응력 텐서 계산
stress_tensor_with_cracks = compute_stress_tensor_with_cracks(
    float(moment_at_x), float(shear_force_at_x),
    width, height, y_val, cracks, x_val
)

# 결과 출력
if stress_tensor_with_cracks is not None:
    print("\n크랙을 고려한 응력 텐서(kPa):")
    for row in stress_tensor_with_cracks:
        formatted_row = [f"{float(value):.2f}" for value in row]
        print("[", " ".join(formatted_row), "]")
    
    # 주응력 계산
    sigma1, sigma2, sigma3 = calculate_principal_stresses_2d_from_3x3(stress_tensor_with_cracks)
    print(f"\n주응력 σ1(kPa): {sigma1:.2f}")
    print(f"주응력 σ2(kPa): {sigma2:.2f}")
    print(f"주응력 σ3(kPa): {sigma3:.2f}")

def visualize_beam_and_sections(bridge_length, height, width, cracks):
    """
    보 전체와 크랙 위치를 시각화하고, 각 크랙의 단면을 별도로 시각화합니다.
    
    Parameters:
        bridge_length (float): 보의 전체 길이 (m)
        height (float): 보의 단면 높이 (m)
        width (float): 보의 단면 폭 (m)
        cracks (list of dict): 크랙 정보 목록, 각 항목은 {'position': x, 'depth': d} 형식
    """
    # 1. 보 전체와 크랙 위치 시각화
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.set_xlim(0, bridge_length)
    ax.set_ylim(-height * 1.5, height * 1.5)

    # 보 그리기
    rect = plt.Rectangle((0, -height / 2), bridge_length, height, color='lightgray', alpha=0.8)
    ax.add_patch(rect)
    ax.axhline(0, color='black', linestyle='--', linewidth=0.8)  # 중립축 표시

    # 크랙 표시
    for crack in cracks:
        crack_pos = crack['x']
        crack_depth = crack['d']
        crack_start_y = -height / 2  # 아랫면에서 시작
        crack_end_y = crack_start_y + crack_depth  # 위로 진행
        ax.plot([crack_pos, crack_pos], [crack_start_y, crack_end_y], color='red', linewidth=2)
        ax.text(crack_pos, crack_end_y + 0.1, f"{crack_depth:.2f}m", color='red', ha='center', fontsize=10)

    # 보 시각화 제목과 라벨
    ax.set_title("Beam and Crack Visualization", fontsize=14)
    ax.set_xlabel("Beam Length (m)", fontsize=12)
    ax.set_ylabel("Height (m)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_aspect('equal', adjustable='box')

    plt.show()

    # 2. 각 크랙의 단면 시각화
    for i, crack in enumerate(cracks, start=1):
        crack_depth = crack['d']
        fig, ax = plt.subplots(figsize=(6, 8))
        
        # 전체 단면 그리기
        rect_section = plt.Rectangle((-width / 2, -height / 2), width, height, color='lightblue', alpha=0.8)
        ax.add_patch(rect_section)

        crack_section = plt.Rectangle((-width / 2, -height / 2), width, crack_depth, color='pink', alpha=0.8)
        ax.add_patch(crack_section)
        
        # 중립축 표시
        ax.axhline(0, color='black', linestyle='--', linewidth=0.8)
        
        # 크랙 표시 (아랫면에서 시작하여 위로 올라감)
        crack_start_y = -height / 2
        crack_end_y = crack_start_y + crack_depth
        ax.plot([0, 0], [crack_start_y, crack_end_y], color='red', linewidth=3)
        ax.text(0.1, crack_end_y, f"{crack_depth:.2f}m", color='red', ha='left', fontsize=10)
        
        # 플롯 설정
        ax.set_title(f"Cross-Section with Crack {i}", fontsize=14)
        ax.set_xlim(-width * 0.6, width * 0.6)
        ax.set_ylim(-height * 0.6, height * 1.2)
        ax.set_xlabel("Width (m)", fontsize=12)
        ax.set_ylabel("Height (m)", fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.set_aspect('equal', adjustable='box')

        plt.tight_layout()
        plt.show()

# 보와 크랙 시각화 실행
visualize_beam_and_sections(bridge_length, height, width, cracks)
#퀘스트2 최종

import numpy as np
import matplotlib.pyplot as plt

def input_material_properties():
    """
    사용자에게 재료 선택과 강도 입력을 요청합니다.
    """
    print("재료를 선택하세요:")
    print("1. 콘크리트")
    print("2. 강철")
    print("3. 목재")
    material_choice = input("선택 (1, 2, 3): ")

    # 재료에 따른 기본 값 및 사용자 입력 요청
    if material_choice == "1":
        print("\n--- 콘크리트 ---")
        tension_strength = float(input("인장 강도 σ_t (kPa)를 입력하세요 (예: 3500): "))
        compression_strength = float(input("압축 강도 σ_c (kPa)를 입력하세요 (예: 30000): "))
        material = "콘크리트"
    elif material_choice == "2":
        print("\n--- 강철 ---")
        tension_strength = float(input("인장 강도 σ_t (kPa)를 입력하세요 (예: 400000): "))
        compression_strength = float(input("압축 강도 σ_c (kPa)를 입력하세요 (예: 400000): "))
        material = "강철"
    elif material_choice == "3":
        print("\n--- 목재 ---")
        tension_strength = float(input("인장 강도 σ_t (kPa)를 입력하세요 (예: 7000): "))
        compression_strength = float(input("압축 강도 σ_c (kPa)를 입력하세요 (예: 30000): "))
        material = "목재"
    else:
        print("잘못된 선택입니다. 프로그램을 종료합니다.")
        return None, None, None

    # 사용자 입력 값 출력
    print("\n--- 선택한 재료 정보 ---")
    print(f"재료: {material}")
    print(f"입력된 인장 강도 σ_t: {tension_strength} kPa")
    print(f"입력된 압축 강도 σ_c: {compression_strength} kPa")

    return tension_strength, compression_strength, material

def evaluate_failure_mohrs_criterion(sigma1, sigma2, sigma_t, sigma_c):
    """
    Mohr's Failure Criterion에 따라 주어진 주응력이 실패 영역에 속하는지 평가합니다.

    Parameters:
        sigma1 (float): 주응력 σ1
        sigma2 (float): 주응력 σ2
        sigma_t (float): 재료의 극한 인장 응력
        sigma_c (float): 재료의 극한 압축 응력

    Returns:
        bool: True이면 실패, False이면 안전
    """
    if sigma1 > 0 and sigma2 > 0:  # 둘 다 인장
        return sigma1 >= sigma_t or sigma2 >= sigma_t
    elif sigma1 < 0 and sigma2 < 0:  # 둘 다 압축
        return sigma1 <= -sigma_c or sigma2 <= -sigma_c
    elif sigma1 > 0 and sigma2 < 0:  # σ1 인장, σ2 압축
        return (sigma1 / sigma_t) + (sigma2 / -sigma_c) >= 1
    elif sigma1 < 0 and sigma2 > 0:  # σ1 압축, σ2 인장
        return (sigma1 / -sigma_c) + (sigma2 / sigma_t) >= 1
    return False

def visualize_mohrs_criterion(sigma1, sigma2, sigma_t, sigma_c):
    """
    Mohr's Failure Criterion을 시각화하고 주어진 주응력을 표시합니다.
    """
    margin = 0.2  # 여유 공간 비율
    sigma_range = np.linspace(-sigma_c * (1 + margin), sigma_t * (1 + margin), 500)
    X, Y = np.meshgrid(sigma_range, sigma_range)

    # Mohr's Criterion 안전 영역 계산
    criterion_safe = np.logical_or.reduce([
        (X > 0) & (Y > 0) & (X < sigma_t) & (Y < sigma_t),  # 둘 다 인장
        (X < 0) & (Y < 0) & (X > -sigma_c) & (Y > -sigma_c),  # 둘 다 압축
        (X > 0) & (Y < 0) & ((X / sigma_t) + (Y / -sigma_c) < 1),  # σ1 인장, σ2 압축
        (X < 0) & (Y > 0) & ((X / -sigma_c) + (Y / sigma_t) < 1)  # σ1 압축, σ2 인장
    ])

    # 초록색 안전 영역(SAFE ZONE) 시각화
    plt.figure(figsize=(8, 8))
    plt.contourf(X, Y, criterion_safe, levels=[0.5, 1], colors=['#c7e9c0'], alpha=0.7)

    # 검은색 경계선 추가 (FAILURE BOUNDARY)
    plt.contour(X, Y, criterion_safe, levels=[0.5], colors='black', linestyles='-', linewidths=1.5)

    # 축 및 제목 설정
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8)  # σ2 축
    plt.axvline(0, color='black', linestyle='--', linewidth=0.8)  # σ1 축
    plt.xlabel('σ1 (kPa)', fontsize=12)
    plt.ylabel('σ2 (kPa)', fontsize=12)
    plt.title("Mohr's Failure Criterion", fontsize=14)

    # 주응력 위치 표시
    plt.scatter(sigma1, sigma2, color='red', s=100, label=f"Point: (σ1={sigma1:.2f}, σ2={sigma2:.2f})")
    plt.legend(loc='upper left', fontsize=10)

    # 주응력 상태 표시
    if evaluate_failure_mohrs_criterion(sigma1, sigma2, sigma_t, sigma_c):
        plt.text(sigma1, sigma2 + 0.5, "FAIL", color='red', fontsize=12, fontweight='bold', ha='center')
    else:
        plt.text(sigma1, sigma2 + 0.5, "SAFE", color='green', fontsize=12, fontweight='bold', ha='center')

    # 그리드와 레이아웃 설정
    plt.grid(linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

# 메인 실행
if __name__ == "__main__":
    tension, compression, material = input_material_properties()
    if tension and compression:
        sigma1 = float(input("\n주응력 σ1 (kPa)를 입력하세요: "))
        sigma2 = float(input("주응력 σ2 (kPa)를 입력하세요: "))

        # Mohr Criterion 평가 및 시각화
        is_failure = evaluate_failure_mohrs_criterion(sigma1, sigma2, tension, compression)
        print(f"\n재료: {material}")
        print(f"Failure Status: {'FAILURE' if is_failure else 'SAFE'}")
        visualize_mohrs_criterion(sigma1, sigma2, tension, compression)

# 프로그램 종료 방지
input("Press Enter to exit...")
