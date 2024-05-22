#include <iostream>
#include <vector>

using namespace std;

typedef vector<vector<int>> Matriz;

void imprimirMatriz(const Matriz &matriz) {
  int n = matriz.size();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < matriz[i].size(); ++j) {
      std::cout << matriz[i][j] << " ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

Matriz sumar(const Matriz &A, const Matriz &B) {
  int n = A.size();
  Matriz C(n, vector<int>(n));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      C[i][j] = A[i][j] + B[i][j];
  return C;
}

Matriz restar(const Matriz &A, const Matriz &B) {
  int n = A.size();
  Matriz C(n, vector<int>(n));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      C[i][j] = A[i][j] - B[i][j];
  return C;
}

Matriz redimensionarMatriz(const Matriz &A, int nuevoTam) {
  int n = A.size();
  Matriz B(nuevoTam, vector<int>(nuevoTam, 0));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      B[i][j] = A[i][j];
  return B;
}

int siguientePotenciaDeDos(int n) {
  int cont = 0;
  if (n && !(n & (n - 1)))
    return n;
  while (n != 0) {
    n >>= 1;
    cont += 1;
  }
  return 1 << cont;
}

Matriz strassen(const Matriz &A, const Matriz &B) {
  int n = A.size();
  if (n == 1) {
    Matriz C(1, vector<int>(1));
    C[0][0] = A[0][0] * B[0][0];
    return C;
  }

  int nuevoTam = n / 2;
  Matriz A11(nuevoTam, vector<int>(nuevoTam));
  Matriz A12(nuevoTam, vector<int>(nuevoTam));
  Matriz A21(nuevoTam, vector<int>(nuevoTam));
  Matriz A22(nuevoTam, vector<int>(nuevoTam));
  Matriz B11(nuevoTam, vector<int>(nuevoTam));
  Matriz B12(nuevoTam, vector<int>(nuevoTam));
  Matriz B21(nuevoTam, vector<int>(nuevoTam));
  Matriz B22(nuevoTam, vector<int>(nuevoTam));

  for (int i = 0; i < nuevoTam; ++i) {
    for (int j = 0; j < nuevoTam; ++j) {
      A11[i][j] = A[i][j];
      A12[i][j] = A[i][j + nuevoTam];
      A21[i][j] = A[i + nuevoTam][j];
      A22[i][j] = A[i + nuevoTam][j + nuevoTam];
      B11[i][j] = B[i][j];
      B12[i][j] = B[i][j + nuevoTam];
      B21[i][j] = B[i + nuevoTam][j];
      B22[i][j] = B[i + nuevoTam][j + nuevoTam];
    }
  }

  Matriz P = strassen(sumar(A11, A22), sumar(B11, B22));
  Matriz Q = strassen(sumar(A21, A22), B11);
  Matriz R = strassen(A11, restar(B12, B22));
  Matriz S = strassen(A22, restar(B21, B11));
  Matriz T = strassen(sumar(A11, A12), B22);
  Matriz U = strassen(restar(A21, A11), sumar(B11, B12));
  Matriz V = strassen(restar(A12, A22), sumar(B21, B22));

  Matriz C(n, vector<int>(n));
  for (int i = 0; i < nuevoTam; ++i) {
    for (int j = 0; j < nuevoTam; ++j) {
      C[i][j] = P[i][j] + S[i][j] - T[i][j] + V[i][j];
      C[i][j + nuevoTam] = R[i][j] + T[i][j];
      C[i + nuevoTam][j] = Q[i][j] + S[i][j];
      C[i + nuevoTam][j + nuevoTam] = P[i][j] + R[i][j] - Q[i][j] + U[i][j];
    }
  }
  return C;
}



Matriz multiplicarStrassen(const Matriz &A, const Matriz &B) {
  int n = A.size();
  int nuevoTam = siguientePotenciaDeDos(n);
  Matriz ARedim = redimensionarMatriz(A, nuevoTam);
  Matriz BRedim = redimensionarMatriz(B, nuevoTam);
  Matriz CRedim = strassen(ARedim, BRedim);
  Matriz C(n, vector<int>(n));

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      C[i][j] = CRedim[i][j];

  return C;
}

int main() {
  int n;
  cout << "Ingresa el tamaño n de las matrices: ";
  cin >> n;

  Matriz A(n, vector<int>(n));
  Matriz B(n, vector<int>(n));

  cout << "Ingresa matriz A:" << endl;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      cin >> A[i][j];

  cout << "Ingresa matriz B:" << endl;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      cin >> B[i][j];

  Matriz C = multiplicarStrassen(A, B);

  cout << "RESULTADO:" << '\n';
  imprimirMatriz(C);

  return 0;
}