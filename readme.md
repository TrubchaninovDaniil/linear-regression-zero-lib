# Linear Regressin (Zero-Dependency)

A high-performance, *zero-dependency* C++ engine for solving linear regression problems. All mathematical operation and decomposition algorithms are implemented from scratch in a dedicated `Maths` class.

The goal of this algorithm is to find the best-fit linear relationship between a set of input features $X$ and targets values $y$. It computes the optimal weight vector $w$ that minimizes the sum of squared residuals, allowing for accurate prediction on new unseen data.

## MSE
We use **Mean Squared Error (MSE)** to measure the quality of the model's predictictions:

$$\large \displaystyle \text{MSE} = \frac{1}{n} \| Y - \hat{Y} \|_2^2$$

To obtain the weights ($w$) that yield the lowest MSE value, we minimize the function by setting its gradient to zero. This leads to the **Normal Equation**:

$$\large \displaystyle (X^T X) w = X^T y$$

## PLU Decomposition
To solve the system $Aw = b$ with maximum numerical stability, we use **LU decomposition with partial pivoting**:

$$\large \displaystyle A = X^T X = PLU$$

$$\text {and}$$

$$\large \displaystyle b = X^T y$$

**where:**
* $P$ - Permutation matrix.
* $L$ - Lower triangular matrix.
* $U$ - Upper triangular matrix.

### **Forward and Back Substitution**
Substituting $A = PLU$ into the system gives $LUw = P^T b$. We solve it as follows: 

We can now solve the transformed system $LUw = P^T b$ in two step:

**Forward Substitution**
We solve for the intermediate vector $z$:
$$\large \displaystyle Lz = P^T b$$

**Backward Substitution**
Using vector $z$, we solve for the final weight vector $w$:
$$\large \displaystyle Uw = z$$

## Quick Start

```bash

# Build the project
g++ -O3 main.cpp maths.cpp -o linreg

# Linux / macOS
./linreg

# Windows
linreg.exe
```

## License
This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.