import sympy as sp

def quaternion_to_rotation_matrix(q):
    """Convert quaternion to rotation matrix symbolically."""
    qw, qx, qy, qz = q
    R = sp.Matrix([
        [1 - 2 * (qy**2 + qz**2), 2 * (qx * qy - qw * qz), 2 * (qx * qz + qw * qy)],
        [2 * (qx * qy + qw * qz), 1 - 2 * (qx**2 + qz**2), 2 * (qy * qz - qw * qx)],
        [2 * (qx * qz - qw * qy), 2 * (qy * qz + qw * qx), 1 - 2 * (qx**2 + qy**2)]
    ])
    return R

def compute_rotation_matrix_pitch(theta):
    """Compute rotation matrix for pitch (y-axis) rotation."""
    R = sp.Matrix([
        [sp.cos(theta), 0, sp.sin(theta)],
        [0, 1, 0],
        [-sp.sin(theta), 0, sp.cos(theta)]
    ])
    return R

def main():
    # Define symbolic parameters
    t = sp.symbols('t')
    l1, lg1, lg2 = sp.symbols('l1 lg1 lg2')  # Link lengths and center of mass distances
    m1, m2 = sp.symbols('m1 m2')  # Mass of links
    g = sp.symbols('g')  # Gravitational acceleration
    I1, Iyy2 = sp.symbols('I1 Iyy2')  # Moments of inertia

    # Generalized coordinates
    qw1, qx1, qy1, qz1 = [sp.Function(f'q1{i}')(t) for i in range(4)]  # Shoulder quaternion
    theta2 = sp.Function('theta2')(t)  # Elbow pitch angle
    theta2_dot = sp.diff(theta2, t)

    # Compute rotation matrices
    R1 = quaternion_to_rotation_matrix([qw1, qx1, qy1, qz1])  # Shoulder rotation
    R2 = compute_rotation_matrix_pitch(theta2)  # Elbow rotation

    # Compute center of mass positions
    O1 = sp.Matrix([0, 0, 0])  # Origin of link 1
    C1 = O1 + R1 * sp.Matrix([lg1, 0, 0])  # Center of mass for link 1
    O2 = O1 + R1 * sp.Matrix([l1, 0, 0])  # Origin of link 2 (end of link 1)
    C2 = O2 + R2 * sp.Matrix([lg2, 0, 0])  # Center of mass for link 2

    # Compute translational kinetic energy
    C1_dot = sp.diff(C1, t)
    C2_dot = sp.diff(C2, t)
# 単位四元数を仮定して、ノルム計算を省略
    T_trans = sp.Rational(1, 2) * m1 * (C1_dot.T * C1_dot)[0] + sp.Rational(1, 2) * m2 * (C2_dot.T * C2_dot)[0]

    # Compute rotational kinetic energy
    # Shoulder quaternion angular velocity
    q1 = sp.Matrix([qw1, qx1, qy1, qz1])
    q1_dot = sp.diff(q1, t)
    omega1 = 2 * sp.Matrix([
        qw1 * q1_dot[1] - qx1 * q1_dot[0] + qy1 * q1_dot[3] - qz1 * q1_dot[2],
        qw1 * q1_dot[2] - qy1 * q1_dot[0] + qz1 * q1_dot[1] - qx1 * q1_dot[3],
        qw1 * q1_dot[3] - qz1 * q1_dot[0] + qx1 * q1_dot[2] - qy1 * q1_dot[1]
    ])
    # Rotational Kinetic Energy の修正
    omega1_norm_squared = omega1.T * omega1  # ノルム二乗計算
    T_rot1 = sp.Rational(1, 2) * I1 * omega1_norm_squared[0]

    T_rot2 = sp.Rational(1, 2) * Iyy2 * theta2_dot**2
    T_rot = T_rot1 + T_rot2

    # Compute potential energy
    U = m1 * g * C1[2] + m2 * g * C2[2]

    # 修正する表記の変換関数
    def simplify_derivatives(expr):
        t = sp.symbols('t')
        functions = [f'q1{i}' for i in range(4)] + ['theta2']  # 対象の変数名
        replacements = {sp.Derivative(sp.Function(func)(t), t): sp.symbols(f'{func}_dot') for func in functions}
        return expr.subs(replacements)

    # エネルギー式の変換後
    U_simplified = simplify_derivatives(U)
    T_trans_simplified = simplify_derivatives(T_trans)
    T_rot_simplified = simplify_derivatives(T_rot)

    # 結果を保存
    with open("energies_simplified.txt", "w") as file:
        file.write("Potential Energy:\n")
        file.write(str(U_simplified) + "\n\n")
        file.write("Translational Kinetic Energy:\n")
        file.write(str(T_trans_simplified) + "\n\n")
        file.write("Rotational Kinetic Energy:\n")
        file.write(str(T_rot_simplified) + "\n")

if __name__ == "__main__":
    main()
