# 20250975 지용건

# 정방행렬을 입력받는 함수
def input_matrix():
    # 정방행렬의 차수 n 입력
    n = int(input("정방행렬의 차수를 입력하세요: "))

    # n이 0 이하이면 잘못된 입력
    if n <= 0:
        print("차수는 양의 정수여야 합니다.")
        return

    print(f"{n}x{n} 정방행렬 A를 입력하세요.")

    # 입력받은 행렬을 저장할 2차원 리스트
    matrix = []

    # n개의 행을 한 행씩 입력
    for i in range(n):
        while True:
            try:
                row_input = input(f"{i+1}행: ").strip()

                # 공백을 기준으로 숫자를 나누어 정수 리스트로 변환
                row_values = [int(x) for x in row_input.split()]

                # 한 행에 정확히 n개의 값이 입력되었는지 확인
                if len(row_values) != n:
                    print(f"입력 오류: 정확히 {n}개의 값을 입력해야 합니다.")
                    continue

                # 입력받은 행을 행렬에 추가
                matrix.append(row_values)
                break

            # 숫자가 아닌 값이 입력된 경우 예외 처리
            except ValueError:
                print("입력 오류: 숫자만 입력해주세요.")

    # 완성된 2차원 리스트 반환
    return matrix


# 행렬식을 계산하는 함수
def determinant(matrix):
    n = len(matrix)

    # 1x1 행렬
    if n == 1:
        return matrix[0][0]

    # 2x2 행렬
    if n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

    det = 0

    # 첫 번째 행을 기준으로 여인수 전개
    for col in range(n):
        minor = []

        for row in matrix[1:]:
            minor.append(row[:col] + row[col+1:])

        sign = (-1) ** col
        det += sign * matrix[0][col] * determinant(minor)

    return det


# 행렬식을 이용하여 역행렬을 계산하는 함수
def inverse_by_determinant(matrix):
    n = len(matrix)
    det = determinant(matrix)

    print("\n행렬식:", det)

    # 행렬식이 0이면 역행렬이 존재하지 않음
    if det == 0:
        print("행렬식이 0이므로 역행렬이 존재하지 않습니다.")
        return

    # 1x1 행렬
    if n == 1:
        return [[1 / matrix[0][0]]]

    # 여인수 행렬 계산
    cofactor_matrix = []

    for i in range(n):
        cofactor_row = []

        for j in range(n):
            minor = []

            # i행과 j열을 제거하여 소행렬 생성
            for row in range(n):
                if row == i:
                    continue

                minor_row = []

                for col in range(n):
                    if col == j:
                        continue

                    minor_row.append(matrix[row][col])

                minor.append(minor_row)

            # 여인수 계산
            cofactor = ((-1) ** (i + j)) * determinant(minor)
            cofactor_row.append(cofactor)

        cofactor_matrix.append(cofactor_row)

    # 여인수 행렬을 전치하여 수반행렬 생성
    adjugate_matrix = []

    for i in range(n):
        row = []

        for j in range(n):
            row.append(cofactor_matrix[j][i])

        adjugate_matrix.append(row)

    # A^-1 = (1 / det(A)) * adj(A)
    inverse_matrix = []

    for i in range(n):
        row = []

        for j in range(n):
            row.append(adjugate_matrix[i][j] / det)

        inverse_matrix.append(row)

    return inverse_matrix


# 가우스-조던 소거법을 이용하여 역행렬을 계산하는 함수
def inverse_by_gauss_jordan(matrix):
    n = len(matrix)

    # 오른쪽에 단위행렬 I를 붙여 [A | I] 형태의 확대행렬 생성
    augmented = []

    for i in range(n):
        row = []

        # 행렬 A 부분
        for j in range(n):
            row.append(float(matrix[i][j]))

        # 단위행렬 I 부분
        for j in range(n):
            if i == j:
                row.append(1.0)
            else:
                row.append(0.0)

        augmented.append(row)

    # 각 열을 차례대로 피벗으로 사용
    for i in range(n):

        # 현재 피벗이 0에 가까우면 아래 행에서
        # 0이 아닌 값을 가진 행을 찾음
        if abs(augmented[i][i]) < 1e-10:
            pivot_row = -1

            for j in range(i + 1, n):
                if abs(augmented[j][i]) >= 1e-10:
                    pivot_row = j
                    break

            # 0이 아닌 피벗을 찾지 못하면 역행렬이 존재하지 않음
            if pivot_row == -1:
                print("가우스-조던 소거법: 역행렬이 존재하지 않습니다.")
                return

            # 현재 행과 피벗으로 사용할 행을 교환
            augmented[i], augmented[pivot_row] = (
                augmented[pivot_row],
                augmented[i]
            )

        # 피벗을 1로 만들기
        pivot = augmented[i][i]

        for j in range(2 * n):
            augmented[i][j] /= pivot

        # 피벗을 제외한 같은 열의 나머지 원소를 모두 0으로 만들기
        for k in range(n):
            if k != i:
                factor = augmented[k][i]

                for j in range(2 * n):
                    augmented[k][j] -= factor * augmented[i][j]

    # [I | A^-1]에서 오른쪽 부분만 추출
    inverse_matrix = []

    for i in range(n):
        row = []

        for j in range(n, 2 * n):
            row.append(round(augmented[i][j], 10))

        inverse_matrix.append(row)

    return inverse_matrix


# 두 방법으로 계산한 역행렬을 비교하는 함수
def compare_matrices(matrix1, matrix2):

    # 둘 중 하나라도 역행렬이 존재하지 않는 경우
    if matrix1 is None or matrix2 is None:
        print("\n역행렬이 존재하지 않아 결과를 비교할 수 없습니다.")
        return

    n = len(matrix1)

    # 두 행렬의 각 원소를 비교
    for i in range(n):
        for j in range(n):

            # 부동소수점 오차를 고려하여 비교
            if abs(matrix1[i][j] - matrix2[i][j]) >= 1e-10:
                print("\n두 방법의 결과가 일치하지 않습니다.")
                return

    print("\n두 방법의 결과가 동일합니다.")


# ---------------- 추가 기능 ----------------

# 두 행렬을 곱하는 함수
def multiply_matrices(matrix1, matrix2):
    n = len(matrix1)
    result = []

    for i in range(n):
        row = []

        for j in range(n):
            value = 0

            # 행렬 곱셈
            for k in range(n):
                value += matrix1[i][k] * matrix2[k][j]

            # 계산값은 그대로 저장
            row.append(value)

        result.append(row)

    return result


# 계산된 역행렬이 올바른지 검산하는 함수
def verify_inverse(matrix, inverse_matrix):

    # 역행렬이 존재하지 않는 경우 검산하지 않음
    if inverse_matrix is None:
        return

    # A × A^-1 계산
    result = multiply_matrices(matrix, inverse_matrix)

    print("\n[추가기능: 역행렬 검산]")
    print("A × A^-1 =")

    # 화면에 출력할 때 부동소수점 오차를 정리
    for row in result:
        print([
            int(round(x))
            if abs(x - round(x)) < 1e-9
            else round(x, 10)
            for x in row
        ])

    n = len(result)

    # 계산 결과가 단위행렬인지 확인
    for i in range(n):
        for j in range(n):

            # 주대각 원소는 1, 나머지는 0이어야 함
            if i == j:
                expected = 1
            else:
                expected = 0

            # 부동소수점 오차를 고려하여 단위행렬과 비교
            if abs(result[i][j] - expected) > 1e-9:
                print("검산 결과: 단위행렬과 일치하지 않습니다.")
                return

    print("검산 결과: 단위행렬과 일치합니다.")
    print("따라서 계산된 역행렬이 올바릅니다.")


# ---------------- main ----------------

def main():

    # 1. 정방행렬 입력
    matrix_a = input_matrix()

    if matrix_a is None:
        return

    print("\n입력된 행렬:")

    for row in matrix_a:
        print(row)

    # 2. 행렬식을 이용한 역행렬 계산
    inverse1 = inverse_by_determinant(matrix_a)

    if inverse1 is not None:
        print("\n[행렬식을 이용한 역행렬]")

        for row in inverse1:
            print([
                int(x) if x.is_integer() else x
                for x in row
            ])

    # 3. 가우스-조던 소거법을 이용한 역행렬 계산
    inverse2 = inverse_by_gauss_jordan(matrix_a)

    if inverse2 is not None:
        print("\n[가우스-조던 소거법을 이용한 역행렬]")

        for row in inverse2:
            print([
                int(x) if x.is_integer() else x
                for x in row
            ])

    # 4. 두 방법으로 계산한 결과 비교
    compare_matrices(inverse1, inverse2)

    # 5. 추가기능 - 계산된 역행렬 검산
    verify_inverse(matrix_a, inverse2)


# 프로그램 실행
if __name__ == "__main__":
    main()