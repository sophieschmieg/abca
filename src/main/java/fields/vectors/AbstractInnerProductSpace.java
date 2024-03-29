package fields.vectors;

import java.lang.reflect.Array;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractModule;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.InnerProductSpace;
import fields.interfaces.Ring;
import fields.interfaces.SesquilinearForm;
import fields.interfaces.ValueField;

public abstract class AbstractInnerProductSpace<T extends Element<T>, S extends Element<S>> extends AbstractModule<T, S>
		implements InnerProductSpace<T, S> {

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public DualVectorSpace<T, S> getDual() {
		return new DualVectorSpace<>(this);
	}

	@Override
	public Ideal<T> annihilator() {
		return getRing().getZeroIdeal();
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	public Ring<T> getRing() {
		return getValueField();
	}

	public Field<T> getField() {
		return getValueField();
	}

	public List<S> getGenerators() {
		return getBasis();
	}

	@Override
	public boolean isLinearIndependent(List<S> s) {
		if (s.size() > dimension()) {
			return false;
		}
		List<Vector<T>> asVectors = new ArrayList<>();
		for (S t : s) {
			asVectors.add(asVector(t));
		}
		Matrix<T> asMatrix = Matrix.fromColumns(asVectors);
		if (getValueField().exactness().equals(Exactness.EXACT)) {
			return asMatrix.getModule(getValueField()).ldup(asMatrix).getSlips().size() == 0;
		}
		return singularValueDecomposition(asMatrix).getRank() == s.size();
	}

	@Override
	public Real valueNorm(S s) {
		Reals r = getValueField().getReals();
		if (s.equals(zero())) {
			return r.zero();
		}
		return r.inverse(inverseValueNorm(s));
	}

	@Override
	public Real inverseValueNorm(S s) {
		Reals r = getValueField().getReals();
		if (dimension() == 1) {
			return r.inverse(getValueField().value(asVector(s).get(1)));
		}
		return r.inverseSqrt(getValueField().value(innerProduct(s, s)));
	}

	public List<S> gramSchmidt(List<S> s) {
		List<S> orthogonal = new ArrayList<>();
		ValueField<T> f = getValueField();
		for (int i = 0; i < s.size(); i++) {
			S vector = s.get(i);
			for (int j = 0; j < i; j++) {
				vector = subtract(vector, scalarMultiply(
						f.divide(innerProduct(s.get(i), s.get(j)), innerProduct(s.get(j), s.get(j))), s.get(j)));
			}
			orthogonal.add(vector);
		}
		return orthogonal;

	}

	public List<S> normedGramSchmidt(List<S> s) {
		List<S> orthogonal = gramSchmidt(s);
		List<S> orthonormal = new ArrayList<>();
		for (S vector : orthogonal) {
			orthonormal.add(normedVector(vector));
		}
		return orthonormal;
	}

	@Override
	public S normedVector(S s) {
		return scalarMultiply(fromReal(inverseValueNorm(s)), s);
	}

	private class BasisExtensionHelper implements Comparable<BasisExtensionHelper> {
		private int n;
		private Real value;

		private BasisExtensionHelper(int n, Real value) {
			this.n = n;
			this.value = value;
		}

		@Override
		public int compareTo(BasisExtensionHelper o) {
			return value.compareTo(o.value);
		}
	}

	@Override
	public List<S> extendToOrthonormalBasis(List<S> s) {
		s = normedGramSchmidt(s);
		if (s.size() == dimension()) {
			return s;
		}
		ValueField<T> f = getValueField();
		Reals r = f.getReals();
		List<Vector<T>> asVectors = new ArrayList<>();
		for (S t : s) {
			asVectors.add(asVector(t));
		}
		List<BasisExtensionHelper> helper = new ArrayList<>();
		for (int i = 0; i < dimension(); i++) {
			Real max = r.zero();
			for (Vector<T> v : asVectors) {
				Real value = f.value(v.get(i + 1));
				if (value.compareTo(max) > 0) {
					max = value;
				}
			}
			helper.add(new BasisExtensionHelper(i + 1, max));
		}
		Collections.sort(helper);
		List<S> basis = new ArrayList<>();
		basis.addAll(s);
		for (int i = 0; i < dimension() - s.size(); i++) {
			basis.add(getUnitVector(helper.get(i).n));
		}
		return normedGramSchmidt(basis);
	}

	@Override
	public SesquilinearForm<T, S> getSesquilinearForm(Matrix<T> t) {
		if (!isHermitean(t)) {
			throw new ArithmeticException("Not hermitean");
		}
		return new SesquilinearForm<>() {
			@Override
			public T evaluate(S t1, S t2) {
				Vector<T> v2 = asVector(t2);
				S p = fromVector(matrixAlgebra().multiply(t, v2));
				return innerProduct(t1, p);
			}
		};
	}

	@Override
	public SesquilinearFormDiagonalizationResult<T, S> diagonalizeSesquilinearForm(SesquilinearForm<T, S> form) {
		return diagonalizeSesquilinearForm(form, Matrix.fromSesquilinearForm(this, form));
	}

	@Override
	public SesquilinearFormDiagonalizationResult<T, S> diagonalizeSesquilinearForm(Matrix<T> matrix) {
		return diagonalizeSesquilinearForm(getSesquilinearForm(matrix), matrix);
	}

	private SesquilinearFormDiagonalizationResult<T, S> diagonalizeSesquilinearForm(SesquilinearForm<T, S> form,
			Matrix<T> matrix) {
		ValueField<T> f = getValueField();
		Reals r = f.getReals();
		Real eps = r.getPowerOfTwo(-r.precision() + 10);
		OrthogonalSimilarResult<T> schur = schurForm(matrix);
		List<Vector<T>> orthonormal = conjugateTranspose(schur.getUnitaryMatrix()).asColumnList();
		List<S> basisPos = new ArrayList<>();
		List<S> basisNeg = new ArrayList<>();
		List<S> basisZero = new ArrayList<>();
		int positive = 0;
		int negative = 0;
		int zero = 0;
		for (Vector<T> v : orthonormal) {
			S vector = fromVector(v);
			T sqr = form.evaluate(vector, vector);
			Real value = f.value(sqr);
			if (value.compareTo(eps) < 0) {
				zero++;
				basisZero.add(vector);
				continue;
			}
			T multiplier = fromReal(r.inverseSqrt(value));
			T normalized = f.multiply(sqr, multiplier, multiplier);
			if (r.value(r.subtract(asComplexNumber(normalized).realPart(), r.one())).compareTo(eps) < 0) {
				positive++;
				basisPos.add(scalarMultiply(multiplier, vector));
			} else {
				negative++;
				basisNeg.add(scalarMultiply(multiplier, vector));
			}
		}
		List<S> basis = new ArrayList<>();
		basis.addAll(basisPos);
		basis.addAll(basisNeg);
		basis.addAll(basisZero);
		return new SesquilinearFormDiagonalizationResult<>(f, basis, positive, negative, zero);
	}

	@Override
	public Matrix<T> conjugateTranspose(Matrix<T> t) {
		List<List<T>> result = new ArrayList<>();
		for (int j = 0; j < t.columns(); j++) {
			List<T> row = new ArrayList<>();
			for (int i = 0; i < t.rows(); i++) {
				row.add(conjugateScalar(t.entry(i + 1, j + 1)));
			}
			result.add(row);
		}
		return new Matrix<>(result);
	}

	private boolean isUpperTriangular(Matrix<T> t) {
		Reals r = getValueField().getReals();
		Real eps = r.power(r.getInteger(2), -r.precision() + 1);
		for (int i = 1; i < t.rows(); i++) {
			for (int j = 0; j < Math.min(i, t.columns()); j++) {
				if (getValueField().value(t.entry(i + 1, j + 1)).compareTo(eps) >= 0) {
					return false;
				}
			}
		}
		return true;
	}
//
//	private boolean isHessenberg(Matrix<T> t) {
//		Reals r = getValueField().getReals();
//		Real eps = r.power(r.getInteger(2), -r.precision() + 1);
//		for (int i = 1; i < t.rows(); i++) {
//			for (int j = 0; j < Math.min(i - 1, t.columns()); j++) {
//				if (getValueField().value(t.entry(i + 1, j + 1)).compareTo(eps) >= 0) {
//					return false;
//				}
//			}
//		}
//		return true;
//	}

	private boolean isHermitean(Matrix<T> t) {
		Matrix<T> conjugateTranspose = conjugateTranspose(t);
		Reals r = getValueField().getReals();
		Real eps = r.power(r.getInteger(2), -r.precision() + 1);
		for (int i = 0; i < t.rows(); i++) {
			for (int j = 0; j < t.columns(); j++) {
				if (getValueField()
						.value(getValueField().subtract(t.entry(i + 1, j + 1), conjugateTranspose.entry(i + 1, j + 1)))
						.compareTo(eps) >= 0) {
					return false;
				}
			}
		}
		return true;
	}

	@Override
	public QRDecompositionResult<T> qrDecomposition(Matrix<T> t) {
		return qrDecomposition(t, false);
	}

	private QRDecompositionResult<T> qrDecompositionHessenberg(Matrix<T> t) {
		if (t.qrResult == null) {
			ValueField<T> f = getValueField();
			Reals r = f.getReals();
			@SuppressWarnings("unchecked")
			T[][] upper = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), t.rows(), t.columns());
			@SuppressWarnings("unchecked")
			T[][] unitary = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), t.rows(), t.rows());
			for (int i = 0; i < t.rows(); i++) {
				for (int j = 0; j < t.columns(); j++) {
					upper[i][j] = i <= j + 1 ? t.entry(i + 1, j + 1) : f.zero();
				}
				for (int j = 0; j < t.rows(); j++) {
					unitary[i][j] = i == j ? f.one() : f.zero();
				}
			}
			// (c -s* 0)
			// U* = (s c* 0)
			// (0 0 1)

			for (int i = 1; i < t.rows(); i++) {
				T diagonal = upper[i - 1][i - 1];
				T subDiagonal = upper[i][i - 1];
				T result = fromReal(
						r.positiveSqrt(asComplexNumber(f.add(f.multiply(conjugateScalar(diagonal), diagonal),
								f.multiply(conjugateScalar(subDiagonal), subDiagonal))).realPart()));
				if (result.equals(f.zero())) {
					continue;
				}
				T cosine = conjugateScalar(f.divide(diagonal, result));
				T sine = conjugateScalar(f.divide(subDiagonal, result));
				// ( c s 0)
				// U = (-s* c* 0)
				// ( 0 0 1)
				// U.A = R
				//
				// U1.A = R1
				// U2.R1 = R2
				// U2.U1.A = R2
				// A = U1*.U2*.R2
				// (a b) (c* -s)
				// (c d) . (s* c)
				upper[i - 1][i - 1] = result;
				upper[i][i - 1] = f.zero();
				for (int j = i; j < t.columns(); j++) {
					T tmp = upper[i - 1][j];
					upper[i - 1][j] = f.add(f.multiply(cosine, upper[i - 1][j]), f.multiply(sine, upper[i][j]));
					upper[i][j] = f.add(f.multiply(conjugateScalar(f.negative(sine)), tmp),
							f.multiply(conjugateScalar(cosine), upper[i][j]));
				}
				for (int j = 0; j < i; j++) {
					T tmp = unitary[j][i - 1];
					unitary[j][i - 1] = f.multiply(unitary[j][i - 1], conjugateScalar(cosine));
					unitary[j][i] = f.multiply(tmp, f.negative(sine));
				}
				unitary[i][i - 1] = conjugateScalar(sine);
				unitary[i][i] = cosine;
			}
			t.qrResult = new QRDecompositionResult<>(new Matrix<>(unitary), new Matrix<>(upper));
		}
		return t.qrResult;
	}

	// (Id - 2ww*)x = x - 2ww*x = x - 2w(w, x)
	private S applyHouseholderLeft(S reflectionVector, S input) {
		T product = innerProduct(reflectionVector, input);
		return subtract(input, scalarMultiply(2, product, reflectionVector));
	}

	private Matrix<T> multiplyHouseholderLeft(S reflectionVector, Matrix<T> t) {
		List<Vector<T>> result = new ArrayList<>();
		for (int i = 0; i < t.columns(); i++) {
			result.add(asVector(applyHouseholderLeft(reflectionVector, fromVector(t.column(i + 1)))));
		}
		return Matrix.fromColumns(result);
	}

	// x^t(Id - 2ww*) = x^t - 2x^tww* = x - 2(bar(x), w)w^*
	private S applyHouseholderRight(S reflectionVector, S input) {
		T product = innerProduct(conjugateVector(input), reflectionVector);
		return subtract(input, scalarMultiply(2, product, conjugateVector(reflectionVector)));
	}

	private Matrix<T> multiplyHouseholderRight(Matrix<T> t, S reflectionVector) {
		List<Vector<T>> result = new ArrayList<>();
		for (int i = 0; i < t.rows(); i++) {
			result.add(asVector(applyHouseholderRight(reflectionVector, fromVector(t.row(i + 1)))));
		}
		return Matrix.fromRows(result);
	}

	@Override
	public QRDecompositionResult<T> qrDecomposition(Matrix<T> t, boolean hessenberg) {
		if (t.columns() != dimension()) {
			throw new ArithmeticException("dimension mismatch!");
		}
		if (t.qrResult != null) {
			return t.qrResult;
		}
		if (hessenberg) {
			return qrDecompositionHessenberg(t);
		}
		Matrix<T> original = t;
		ValueField<T> f = getValueField();
		Reals r = f.getReals();
		MatrixAlgebra<T> algebra = matrixAlgebra();
		Matrix<T> unitary = algebra.one();
		Complex c = Complex.c(r.precision());
		for (int j = 0; j < t.columns() - 1; j++) {
			T sum = f.zero();
			List<T> columnValues = new ArrayList<>();
			for (int i = 0; i < j; i++) {
				columnValues.add(f.zero());
			}
			for (int i = j; i < t.rows(); i++) {
				columnValues.add(t.entry(i + 1, j + 1));
				sum = f.add(f.multiply(conjugateScalar(t.entry(i + 1, j + 1)), t.entry(i + 1, j + 1)), sum);
			}
			S column = fromVector(new Vector<>(columnValues));
			Real norm = r.positiveSqrt(asComplexNumber(sum).realPart());
			T pivot = t.entry(j + 1, j + 1);
			ComplexNumber complexPivot = asComplexNumber(pivot);
			ComplexNumber pivotValue = c.getEmbedding(c.value(complexPivot));
			ComplexNumber normedPivot = pivotValue.equals(c.zero()) ? c.one() : c.divide(complexPivot, pivotValue);
			T alpha = f.multiply(-1, fromComplexNumber(normedPivot), fromReal(norm));
			S reflectionVector = subtract(column, scalarMultiply(alpha, getUnitVector(j + 1)));
			if (reflectionVector.equals(zero())) {
				continue;
			}
			reflectionVector = normedVector(reflectionVector);
			unitary = multiplyHouseholderRight(unitary, reflectionVector);
			t = multiplyHouseholderLeft(reflectionVector, t);
		}
		original.qrResult = new QRDecompositionResult<>(unitary, forceZero(t));
		return original.qrResult;
	}

	@Override
	public OrthogonalSimilarResult<T> hessenbergForm(Matrix<T> t) {
		if (t.columns() != dimension() || t.rows() != t.columns()) {
			throw new ArithmeticException("dimension mismatch!");
		}
		if (t.rows() <= 2) {
			return new OrthogonalSimilarResult<>(matrixAlgebra().one(), t);
		}
		ValueField<T> f = getValueField();
		Reals r = f.getReals();
		Complex c = Complex.c(r.precision());
		MatrixAlgebra<T> algebra = matrixAlgebra();
		Matrix<T> unitary = algebra.one();
		for (int j = 0; j < t.columns() - 2; j++) {
			T sum = f.zero();
			List<T> columnValues = new ArrayList<>();
			for (int i = 0; i <= j; i++) {
				columnValues.add(f.zero());
			}
			for (int i = j + 1; i < t.rows(); i++) {
				columnValues.add(t.entry(i + 1, j + 1));
				sum = f.add(f.multiply(conjugateScalar(t.entry(i + 1, j + 1)), t.entry(i + 1, j + 1)), sum);
			}
			S column = fromVector(new Vector<>(columnValues));
			Real norm = r.positiveSqrt(c.asReal(asComplexNumber(sum)));
			T pivot = t.entry(j + 2, j + 1);
			ComplexNumber complexPivot = asComplexNumber(pivot);
			ComplexNumber pivotValue = c.getEmbedding(c.value(complexPivot));
			ComplexNumber normedPivot = pivotValue.equals(c.zero()) ? c.getInteger(-1)
					: c.divide(c.conjugate(complexPivot), pivotValue);
			S reflectionVector = add(scalarMultiply(fromComplexNumber(normedPivot), column),
					scalarMultiply(fromReal(norm), getUnitVector(j + 2)));
			if (reflectionVector.equals(zero())) {
				continue;
			}
			reflectionVector = normedVector(reflectionVector);
			t = multiplyHouseholderLeft(reflectionVector, multiplyHouseholderRight(t, reflectionVector));
			unitary = multiplyHouseholderRight(unitary, reflectionVector);
		}
		// w - 2ww*w = w - 2w = -w
//		Matrix<T> hessenberg = algebra.multiply(conjugateTranspose(unitary), t, unitary);
//		if (isHessenberg(hessenberg)) {
//			return new OrthogonalSimilarResult<>(unitary, hessenberg);
//		}
//		InnerProductSpace<T, S> subSpace = withDimension(dimension() - 1);
//		OrthogonalSimilarResult<T> decomposition = subSpace
//				.hessenbergForm(hessenberg.subMatrix(2, t.rows(), 2, t.columns()));
//		List<Matrix<T>> decomposedUnitary = new ArrayList<>();
//		decomposedUnitary.add(new Matrix<>(Collections.singletonList(Collections.singletonList(f.one()))));
//		decomposedUnitary.add(decomposition.getUnitaryMatrix());
//		Matrix<T> lowerUnitary = Matrix.fromBlockDiagonalMatrix(decomposedUnitary, f.zero());
//		Matrix<T> resultUnitary = algebra.multiply(unitary, lowerUnitary);
//		Matrix<T> hessenbergForm = forceZero(algebra.multiply(conjugateTranspose(resultUnitary), t, unitary));
		return new OrthogonalSimilarResult<>(unitary, forceZero(t));
	}

	private Matrix<T> forceZero(Matrix<T> t) {
		ValueField<T> f = getValueField();
		Reals r = f.getReals();
		int precision = dimension();
		for (int i = 0; i < t.rows(); i++) {
			for (int j = 0; j < t.columns(); j++) {
				IntE value = f.value(t.entry(i + 1, j + 1)).roundUp();
				precision = Math.max(precision, value.getValue().bitLength());
			}
		}
		Real eps = r.power(r.getInteger(2), -r.precision() + precision);
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < t.rows(); i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < t.columns(); j++) {
				T entry = t.entry(i + 1, j + 1);
				if (f.value(entry).compareTo(eps) < 0) {
					row.add(getValueField().zero());
				} else {
					ComplexNumber asComplex = asComplexNumber(entry);
					if (r.value(asComplex.complexPart()).compareTo(eps) < 0) {
						row.add(fromReal(asComplex.realPart()));
					} else {
						row.add(entry);
					}

				}
			}
			result.add(row);
		}
		return new Matrix<>(result);
	}

	@Override
	public OrthogonalSimilarResult<T> schurForm(Matrix<T> t) {
		if (t.schurFormResult == null) {
			MatrixAlgebra<T> algebra = matrixAlgebra();
			if (t.rows() == 1 && dimension() == 1) {
				return new OrthogonalSimilarResult<>(algebra.one(), t);
			}
			OrthogonalSimilarResult<T> hessenberg = hessenbergForm(t);
			Matrix<T> it = hessenberg.getOrthogonallySimilarMatrix();
			Matrix<T> unitary = hessenberg.getUnitaryMatrix();
			// it = u^* t u
			ValueField<T> f = getValueField();
			Reals r = f.getReals();
			Complex c = Complex.c(r.precision());
			Real eps = r.power(r.getInteger(2), -r.precision() + dimension());
			int limit = t.rows();
			int counter = 0;
			while (true) {
				counter++;
				if (counter > 1000) {
					throw new ArithmeticException("Schurr Form not converging!");
				}
				// System.out.println(algebra.multiply(algebra.multiply(unitary,
				// it),conjugateTranspose( unitary)));
//				ComplexNumber norm = asComplexNumber(
//						f.subtract(f.multiply(it.entry(limit - 1, limit - 1), it.entry(limit, limit)),
//								f.multiply(it.entry(limit - 1, limit), it.entry(limit, limit - 1))));
				ComplexNumber trace = asComplexNumber(f.add(it.entry(limit - 1, limit - 1), it.entry(limit, limit)));
				ComplexNumber diagonalDiff = asComplexNumber(
						f.divide(f.subtract(it.entry(limit - 1, limit - 1), it.entry(limit, limit)), f.getInteger(2)));
				ComplexNumber counterDiagonalProduct = asComplexNumber(
						f.multiply(it.entry(limit, limit - 1), it.entry(limit - 1, limit)));
				ComplexNumber discriminant = c.add(c.multiply(diagonalDiff, diagonalDiff), counterDiagonalProduct);
				ComplexNumber root;
				if (!c.isReal(discriminant) || discriminant.realPart().compareTo(r.zero()) < 0) {
					root = c.roots(discriminant, 2).keySet().iterator().next();
				} else {
					root = c.getEmbedding(r.positiveSqrt(discriminant.realPart()));
				}
				T option1 = fromComplexNumber(c.add(c.divide(trace, c.getInteger(2)), root));
				T option2 = fromComplexNumber(c.subtract(c.divide(trace, c.getInteger(2)), root));
				Real option1Value = c.norm(asComplexNumber(f.subtract(it.entry(limit, limit), option1)));
				Real option2Value = c.norm(asComplexNumber(f.subtract(it.entry(limit, limit), option2)));
				T shift = option1Value.compareTo(option2Value) > 0 ? option2 : option1;
				Matrix<T> shiftDiagonal = algebra.scalarMultiply(shift, algebra.one());
				QRDecompositionResult<T> qr = qrDecomposition(algebra.subtract(it, shiftDiagonal), true);
//				System.out.println(algebra.subtract(it, shiftDiagonal));
//				System.out.println(algebra.multiply(qr.getUnitaryMatrix(), qr.getUpperTriangularMatrix()));
//				System.out.println(algebra.multiply(conjugateTranspose(qr.getUnitaryMatrix()), qr.getUnitaryMatrix()));
				// q r = it
				// q r = u^* t u
				// it' = r q
				// it' = q^* q r q = q^* it q = q^* u^* t u q
				it = algebra.add(algebra.multiply(qr.getUpperTriangularMatrix(), qr.getUnitaryMatrix()), shiftDiagonal);
				unitary = algebra.multiply(unitary, qr.getUnitaryMatrix());
				if (f.value(it.entry(limit, limit - 1)).compareTo(eps) < 0) {
					it = forceZero(it);
					limit--;
				}
				if (isUpperTriangular(it)) {
					break;
				}
			}
			t.schurFormResult = new OrthogonalSimilarResult<>(conjugateTranspose(unitary), it);
		}
		return t.schurFormResult;
	}

	private class SVDHelper implements Comparable<SVDHelper> {
		private int index;
		private Real value;

		private SVDHelper(int index, Real value) {
			this.index = index;
			this.value = value;
		}

		@Override
		public int compareTo(SVDHelper o) {
			return -value.compareTo(o.value);
		}

	}

	@Override
	public SingularValueDecompositionResult<T> singularValueDecomposition(Matrix<T> t) {
		if (t.svdResult == null) {
			ValueField<T> f = getValueField();
			Reals r = f.getReals();
			InnerProductSpace<T, S> columnSpace = withDimension(t.columns());
			InnerProductSpace<T, S> rowSpace = withDimension(t.rows());
			MatrixModule<T> module = t.getModule(getValueField());
			OrthogonalSimilarResult<T> domainSchurr = columnSpace
					.schurForm(module.domainAlgebra().multiply(conjugateTranspose(t), t));
			List<List<T>> singularValues = new ArrayList<>();
			List<T> pseudoInverseDiagonal = new ArrayList<>();
			Real eps = r.power(r.getInteger(2), -r.precision() / 2);
			int rank = 0;
			List<SVDHelper> helper = new ArrayList<>();
			for (int i = 0; i < t.columns(); i++) {
				helper.add(new SVDHelper(i, f.value(domainSchurr.getOrthogonallySimilarMatrix().entry(i + 1, i + 1))));
			}
			Collections.sort(helper);
			int[] permutation = new int[t.columns()];
			for (int i = 0; i < t.columns(); i++) {
				permutation[i] = helper.get(i).index;
			}
			for (int i = 0; i < t.rows(); i++) {
				List<T> row = new ArrayList<>();
				for (int j = 0; j < t.columns(); j++) {
					Real diagonal = i == j ? helper.get(i).value : r.zero();
					if (i == j && diagonal.compareTo(eps) > 0) {
						rank++;
						T sqrt = fromReal(r.positiveSqrt(diagonal));
						row.add(sqrt);
						pseudoInverseDiagonal.add(f.inverse(sqrt));
					} else {
						row.add(f.zero());
					}
				}
				singularValues.add(row);
			}
			Matrix<T> rightUnitary = module.domainAlgebra().permuteRows(permutation, domainSchurr.getUnitaryMatrix());
			Matrix<T> diagonal = new Matrix<>(singularValues);
			List<S> leftColumns = new ArrayList<>();
			Matrix<T> halfInverted = module.multiply(t, conjugateTranspose(rightUnitary));
			for (int i = 0; i < rank; i++) {
				leftColumns.add(rowSpace.scalarMultiply(pseudoInverseDiagonal.get(i),
						rowSpace.fromVector(halfInverted.column(i + 1))));
			}
			leftColumns = rowSpace.extendToOrthonormalBasis(leftColumns);
			List<Vector<T>> leftColumnVectors = new ArrayList<>();
			for (S column : leftColumns) {
				leftColumnVectors.add(rowSpace.asVector(column));
			}
			Matrix<T> leftUnitary = Matrix.fromColumns(leftColumnVectors);
			t.svdResult = new SingularValueDecompositionResult<>(rank, conjugateTranspose(leftUnitary), diagonal,
					rightUnitary);
		}
		return t.svdResult;
	}

	@Override
	public Matrix<T> pseudoInverse(Matrix<T> t) {
		if (t.pseudoInverse == null) {
			ValueField<T> f = getValueField();
			SingularValueDecompositionResult<T> svd = singularValueDecomposition(t);
			List<List<T>> diagonal = new ArrayList<>();
			for (int i = 0; i < t.columns(); i++) {
				List<T> row = new ArrayList<>();
				for (int j = 0; j < t.rows(); j++) {
					if (i == j && !svd.getDiagonalMatrix().entry(j + 1, i + 1).equals(f.zero())) {
						row.add(f.inverse(svd.getDiagonalMatrix().entry(j + 1, i + 1)));
					} else {
						row.add(f.zero());
					}
				}
				diagonal.add(row);
			}
			Matrix<T> pseudoInverseDiagonal = new Matrix<>(diagonal);
			MatrixModule<T> mm = pseudoInverseDiagonal.getModule(f);
			t.pseudoInverse = mm.multiply(
					mm.multiply(conjugateTranspose(svd.getRightUnitaryMatrix()), pseudoInverseDiagonal),
					svd.getLeftUnitaryMatrix());
		}
		return t.pseudoInverse;
	}

	@Override
	public Optional<Real> conditionNumber(Matrix<T> t) {
		SingularValueDecompositionResult<T> svd = singularValueDecomposition(t);
		Real max = null;
		Real min = null;
		for (int i = 0; i < Math.min(t.rows(), t.columns()); i++) {
			Real sigma = getValueField().value(svd.getDiagonalMatrix().entry(i + 1, i + 1));
			if (max == null || sigma.compareTo(max) > 0) {
				max = sigma;
			}
			if (min == null || sigma.compareTo(min) < 0) {
				min = sigma;
			}
		}
		Real zero = getValueField().getReals().zero();
		if (min.equals(zero)) {
			return Optional.empty();
		}
		return Optional.of(getValueField().getReals().divide(max, min));
	}

	@Override
	public boolean isSubModuleMemberModule(MatrixModule<T> module, Matrix<T> m, Vector<T> b) {
		ValueField<T> f = getValueField();
		int precisionDiscount = 0;
		for (int i = 0; i < m.rows(); i++) {
			for (int j = 0; j < m.columns(); j++) {
				int exp = f.value(m.entry(i + 1, j + 1)).exponent();
				if (precisionDiscount < exp) {
					precisionDiscount = exp;
				}
			}
		}
		for (int i = 0; i < b.dimension(); i++) {
			int exp = f.value(b.get(i + 1)).exponent();
			if (precisionDiscount < exp) {
				precisionDiscount = exp;
			}
		}
		Real eps = f.getReals().getPowerOfTwo(-f.getReals().precision() + precisionDiscount + 10);
		SingularValueDecompositionResult<T> svd = singularValueDecomposition(m);
		Vector<T> x = module.codomainAlgebra().multiply(svd.getLeftUnitaryMatrix(), b);
		for (int i = svd.getRank(); i < m.rows(); i++) {
			if (f.value(x.get(i + 1)).compareTo(eps) >= 0) {
				return false;
			}
		}
		return true;
	}

	@Override
	public Vector<T> asSubModuleMemberModule(MatrixModule<T> module, Matrix<T> m, Vector<T> b) {
		Matrix<T> pseudoInverse = pseudoInverse(m);
		return pseudoInverse.getModule(getField()).multiply(pseudoInverse, b);
	}

	@Override
	public List<Vector<T>> syzygyProblemModule(MatrixModule<T> module, Matrix<T> m) {
		SingularValueDecompositionResult<T> svd = singularValueDecomposition(m);
		Matrix<T> inverted = conjugateTranspose(svd.getRightUnitaryMatrix());
		List<Vector<T>> result = new ArrayList<>();
		for (int i = svd.getRank(); i < m.columns(); i++) {
			result.add(inverted.column(i + 1));
		}
		return result;
	}

	@Override
	public List<Vector<T>> simplifySubModuleGeneratorsModule(MatrixModule<T> module, Matrix<T> m) {
		SingularValueDecompositionResult<T> svd = singularValueDecomposition(m);
		List<Vector<T>> result = new ArrayList<>();
		Matrix<T> left = conjugateTranspose(svd.getLeftUnitaryMatrix());
		for (int i = 0; i < svd.getRank(); i++) {
			result.add(left.column(i + 1));
		}
		return result;
	}
}
