package fields.vectors;

import java.lang.reflect.Array;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.InnerProductSpace;
import fields.interfaces.InnerProductSpace.QRDecompositionResult;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;

public class MatrixAlgebra<T extends Element<T>> extends AbstractAlgebra<T, Matrix<T>> {
	private Ring<T> ring;
	private int dimension;
	private Matrix<T> id;
	private MatrixModule<T> module;

	MatrixAlgebra(Ring<T> ring, int dimension, FreeModule<T> free) {
		this.ring = ring;
		this.dimension = dimension;
		this.module = new MatrixModule<T>(ring, dimension, free, this);
		List<List<T>> id = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			List<T> idRow = new ArrayList<>();
			for (int j = 0; j < dimension; j++) {
				if (i == j) {
					idRow.add(ring.one());
				} else {
					idRow.add(ring.zero());
				}
			}
			id.add(idRow);
		}
		this.id = new Matrix<T>(id);
	}

	@Override
	public Exactness exactness() {
		return ring.exactness();
	}

	@Override
	public Ring<T> getRing() {
		return ring;
	}

	@Override
	public Matrix<T> zero() {
		return module.zero();
	}

	@Override
	public Matrix<T> add(Matrix<T> s1, Matrix<T> s2) {
		return module.add(s1, s2);
	}

	@Override
	public Matrix<T> negative(Matrix<T> s) {
		return module.negative(s);
	}

	@Override
	public Matrix<T> getEmbedding(T t) {
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < dimension; j++) {
				if (i == j) {
					row.add(t);
				} else {
					row.add(ring.zero());
				}
			}
			result.add(row);
		}
		return new Matrix<T>(result);
	}

	@Override
	public Matrix<T> scalarMultiply(T t, Matrix<T> s) {
		return module.scalarMultiply(t, s);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public Ideal<T> annihilator() {
		return ring.getZeroIdeal();
	}
	
	@Override
	public List<Vector<T>> getSyzygies() {
		return Collections.emptyList();
	}

	@Override
	public Matrix<T> getRandomElement() {
		return module.getRandomElement();
	}

	@Override
	public boolean isFinite() {
		return module.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return module.getNumberOfElements();
	}

	@Override
	public Iterator<Matrix<T>> iterator() {
		return module.iterator();
	}

	@Override
	public Matrix<T> one() {
		return id;
	}

	@Override
	public BigInteger characteristic() {
		return ring.characteristic();
	}

	@Override
	public Matrix<T> multiply(Matrix<T> t1, Matrix<T> t2) {
		return module.multiply(t1, t2);
	}

	public Matrix<T> swapOperation(int iIndex, int jIndex) {
		if (iIndex == jIndex) {
			return one();
		}
		iIndex--;
		jIndex--;
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(one().entry(1, 1).getClass(), dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				if ((i == iIndex && j == iIndex) || (i == jIndex && j == jIndex)) {
					matrix[i][j] = ring.zero();
				} else if ((i == iIndex && j == jIndex) || (i == jIndex && j == iIndex)) {
					matrix[i][j] = ring.one();
				} else if (i == j) {
					matrix[i][j] = ring.one();
				} else {
					matrix[i][j] = ring.zero();
				}
			}
		}
		return new Matrix<T>(matrix);
	}

	public Matrix<T> rowShearOperation(int iIndex, int jIndex, T value) {
		return module.rowShearOperation(iIndex, jIndex, value);
	}

	public Matrix<T> colShearOperation(int iIndex, int jIndex, T value) {
		return module.colShearOperation(iIndex, jIndex, value);
	}

	public MatrixModule<T>.LDUPResult ldup(Matrix<T> m) {
		return module.ldup(m);
	}

	public Matrix<T> lupUpperTriangle(Matrix<T> m) {
		return module.lupUpperTriangle(m);
	}

	public MatrixModule<T>.SmithNormalFormResult smithNormalForm(Matrix<T> t) {
		return module.smithNormalForm(t);
	}

	public Matrix<T> adjugate(Matrix<T> t) {
		MatrixModule<T>.SmithNormalFormResult gauss = smithNormalForm(t);
		return multiply(gauss.getColOperationsInverse(), gauss.getRowOperationsInverse());
	}

	public Set<T> eigenValues(Matrix<T> t) {
		if (!(ring instanceof Field<?>)) {
			throw new UnsupportedOperationException("Not over a Vector Space");
		}
		UnivariatePolynomial<T> charPoly = ring.getUnivariatePolynomialRing().normalize(characteristicPolynomial(t));
		return ring.roots(charPoly).keySet();
	}

	public SubVectorSpace<T, Vector<T>> eigenSpace(Matrix<T> t, T eigenValue) {
		if (!(ring instanceof Field<?>)) {
			throw new UnsupportedOperationException("Not over a Vector Space");
		}
		return eigenSpace(t, eigenValue, new FiniteVectorSpace<>((Field<T>) ring, dimension));
	}

	private SubVectorSpace<T, Vector<T>> eigenSpace(Matrix<T> t, T eigenValue, FiniteVectorSpace<T> space) {
		Matrix<T> mod = subtract(t, scalarMultiply(eigenValue, one()));
		List<Vector<T>> eigenSpace = kernelBasis(mod);
		return new SubVectorSpace<>(space, eigenSpace);
	}

	public Map<T, SubVectorSpace<T, Vector<T>>> eigenSpaces(Matrix<T> t) {
		if (!(ring instanceof Field<?>)) {
			throw new UnsupportedOperationException("Not over a Vector Space");
		}
		FiniteVectorSpace<T> space = new FiniteVectorSpace<>((Field<T>) ring, dimension);
		Map<T, SubVectorSpace<T, Vector<T>>> result = new TreeMap<>();
		Set<T> eigenValues = eigenValues(t);
		for (T eigenValue : eigenValues) {
			result.put(eigenValue, eigenSpace(t, eigenValue, space));
		}
		return result;
	}

	public static class DiagonalizationResult<T extends Element<T>> {
		private List<Vector<T>> basis;
		private Matrix<T> diagonalMatrix;
		private Matrix<T> baseChange;
		private Matrix<T> inverseBaseChange;

		private DiagonalizationResult(List<Vector<T>> basis, Matrix<T> diagonalMatrix, Matrix<T> baseChange,
				Matrix<T> inverseBaseChange) {
			this.basis = basis;
			this.diagonalMatrix = diagonalMatrix;
			this.baseChange = baseChange;
			this.inverseBaseChange = inverseBaseChange;
		}

		public List<Vector<T>> getBasis() {
			return basis;
		}

		public Matrix<T> getDiagonalMatrix() {
			return diagonalMatrix;
		}

		public Matrix<T> getBaseChange() {
			return baseChange;
		}

		public Matrix<T> getInverseBaseChange() {
			return inverseBaseChange;
		}
	}

	public Optional<DiagonalizationResult<T>> diagonalize(Matrix<T> t) {
		if (!(ring instanceof Field<?>)) {
			throw new UnsupportedOperationException("Not over a Vector Space");
		}
		Map<T, SubVectorSpace<T, Vector<T>>> eigenSpaces = eigenSpaces(t);
		List<Vector<T>> basis = new ArrayList<>();
		List<T> diagonal = new ArrayList<>();
		for (T eigenValue : eigenSpaces.keySet()) {
			SubVectorSpace<T, Vector<T>> eigenSpace = eigenSpaces.get(eigenValue);
			for (int i = 0; i < eigenSpace.dimension(); i++) {
				diagonal.add(eigenValue);
			}
			basis.addAll(eigenSpace.getBasis());
		}
		if (basis.size() != dimension) {
			return Optional.empty();
		}
		Matrix<T> baseChange = Matrix.fromColumns(basis);
		Matrix<T> inverseBaseChange = inverse(baseChange);
		List<List<T>> diagonalMatrixList = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < dimension; j++) {
				if (i == j) {
					row.add(diagonal.get(i));
				} else {
					row.add(ring.zero());
				}
			}
			diagonalMatrixList.add(row);
		}
		Matrix<T> diagonalMatrix = new Matrix<>(diagonalMatrixList);
		return Optional.of(new DiagonalizationResult<>(basis, diagonalMatrix, baseChange, inverseBaseChange));
	}

	public UnivariatePolynomial<T> characteristicPolynomial(Matrix<T> t) {
		UnivariatePolynomialRing<T> r = ring.getUnivariatePolynomialRing();
		MatrixAlgebra<Polynomial<T>> a = new FreeModule<>(r, dimension).matrixAlgebra();
		List<List<Polynomial<T>>> asPolynomial = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			asPolynomial.add(new ArrayList<>());
			for (int j = 0; j < dimension; j++) {
				asPolynomial.get(i).add(r.getEmbedding(t.entry(i + 1, j + 1)));
			}
		}
		return r.toUnivariate(a.determinant(a.subtract(a.getEmbedding(r.getVar(1)), new Matrix<>(asPolynomial))));
	}

	public T determinant(Matrix<T> t) {
		MatrixModule<T>.LDUPResult ldup = module.ldup(t);
		return ring.multiply(ldup.getSignum(), ldup.getDeterminant());
	}

	public T trace(Matrix<T> t) {
		T tr = ring.zero();
		for (int i = 0; i < dimension; i++) {
			tr = ring.add(tr, t.entry(i + 1, i + 1));
		}
		return tr;
	}

	public Vector<T> multiply(Matrix<T> t, Vector<T> x) {
		return module.multiply(t, x);
	}

	public Vector<T> multiply(int[] permutation, Vector<T> x) {
		return module.multiply(permutation, x);
	}

	public List<Vector<T>> kernelBasis(Matrix<T> t) {
		return module.kernelBasis(t);
	}

	public boolean hasSolution(Matrix<T> t, Vector<T> b) {
		return module.hasSolution(t, b);
	}

	public Vector<T> solve(Matrix<T> t, Vector<T> b) {
		return module.solve(t, b);
	}

	@Override
	public boolean isUnit(Matrix<T> t) {
		return ring.isUnit(determinant(t));
	}

	@SuppressWarnings("unchecked")
	@Override
	public Matrix<T> inverse(Matrix<T> t) {
		if (t.pseudoInverse == null) {
			if (module.hasInnerProductSpace()) {
				InnerProductSpace<T, Vector<T>> space = module.domainInnerProductSpace();
				QRDecompositionResult<T> qr = space.qrDecomposition(t);
				t.pseudoInverse = multiply(invertUpperTriangleMatrix(qr.getUpperTriangularMatrix()),
						space.conjugateTranspose(qr.getUnitaryMatrix()));
			} else {
				MatrixModule<T>.SmithNormalFormResult gauss = smithNormalForm(t);
				T[][] matrix = (T[][]) Array.newInstance(t.entry(1, 1).getClass(), dimension, dimension);
				for (int i = 0; i < dimension; i++) {
					for (int j = 0; j < dimension; j++) {
						if (i == j) {
							matrix[i][j] = ring.inverse(gauss.getDiagonalMatrix().entry(i + 1, j + 1));
						} else {
							matrix[i][j] = ring.zero();
						}
					}
				}
				t.pseudoInverse = multiply(gauss.getColOperationsInverse(), new Matrix<T>(matrix),
						gauss.getRowOperationsInverse());
			}
		}
		return t.pseudoInverse;
	}

	public Matrix<T> invertUpperTriangleMatrix(Matrix<T> t) {
		return module.invertUpperTriangleMatrix(t);
	}

	public Matrix<T> invertLowerTriangleMatrix(Matrix<T> lowerTriangle) {
		return module.invertLowerTriangleMatrix(lowerTriangle);
	}

	public int[] invertPermutation(int[] permutation) {
		return module.invertPermutation(permutation);
	}

	public Matrix<T> permuteRows(int[] permutation, Matrix<T> t) {
		return module.permuteRows(permutation, t);
	}

	public Matrix<T> permuteColumns(int[] permutation, Matrix<T> t) {
		return module.permuteColumns(permutation, t);
	}

	@Override
	public boolean isIntegral() {
		return dimension == 0 && ring.isIntegral();
	}

	@Override
	public boolean isReduced() {
		return dimension == 0 && ring.isReduced();
	}

	@Override
	public boolean isIrreducible() {
		return dimension == 0 && ring.isIrreducible();
	}

	@Override
	public boolean isZeroDivisor(Matrix<T> t) {
		return kernelBasis(t).size() == 0;
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	@Override
	public FactorizationResult<Matrix<T>, Matrix<T>> uniqueFactorization(Matrix<T> t) {
		throw new ArithmeticException("Not a UFD");
	}

	@Override
	public boolean isCommutative() {
		return dimension == 1;
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return false;
	}

	@Override
	public boolean isDedekindDomain() {
		return false;
	}

	@Override
	public boolean isDivisible(Matrix<T> dividend, Matrix<T> divisor) {
		throw new UnsupportedOperationException();
	}

	@Override
	public QuotientAndRemainderResult<Matrix<T>> quotientAndRemainder(Matrix<T> dividend, Matrix<T> divisor) {
		throw new RuntimeException();
	}

	@Override
	public BigInteger euclidMeasure(Matrix<T> t) {
		return null;
	}

	@Override
	public Matrix<T> projectToUnit(Matrix<T> t) {
		return one();
	}

	@Override
	public Iterable<Matrix<T>> getUnits() {
		return new Iterable<Matrix<T>>() {

			@Override
			public Iterator<Matrix<T>> iterator() {
				return new Iterator<Matrix<T>>() {
					private Matrix<T> next = null;
					private Iterator<Matrix<T>> it = MatrixAlgebra.this.iterator();

					private void setNext() {
						if (next == null) {
							while (it.hasNext()) {
								Matrix<T> c = it.next();
								if (MatrixAlgebra.this.isUnit(c)) {
									next = c;
									break;
								}
							}
						}
					}

					@Override
					public boolean hasNext() {
						setNext();
						return next != null;
					}

					@Override
					public Matrix<T> next() {
						setNext();
						Matrix<T> result = next;
						next = null;
						setNext();
						return result;
					}
				};
			}
		};
	}

	public int rank(Matrix<T> m) {
		return module.rank(m);
	}

	public MatrixModule<T> asModule() {
		return module;
	}

	public Matrix<T> transpose(Matrix<T> t) {
		return module.transpose(t);
	}

	@Override
	public List<Matrix<T>> getModuleGenerators() {
		return module.getModuleGenerators();
	}

	@Override
	public Vector<T> asVector(Matrix<T> m) {
		return module.asVector(m);
	}

	public FreeModule<T> asFreeModule() {
		return module.asFreeModule();
	}

	@Override
	public boolean isLinearIndependent(List<Matrix<T>> s) {
		return module.isLinearIndependent(s);
	}

	@Override
	public List<Vector<T>> nonTrivialCombinations(List<Matrix<T>> s) {
		return module.nonTrivialCombinations(s);
	}

	@Override
	public boolean isGeneratingModule(List<Matrix<T>> s) {
		return module.isGeneratingModule(s);
	}

	@Override
	public IdealResult<Matrix<T>, ?> getIdealWithTransforms(List<Matrix<T>> generators) {
		throw new ArithmeticException("Not a commutative ring!");
	}

	@Override
	public Ideal<Matrix<T>> intersect(Ideal<Matrix<T>> t1, Ideal<Matrix<T>> t2) {
		throw new ArithmeticException("Not a commutative ring!");
	}

	@Override
	public Ideal<Matrix<T>> radical(Ideal<Matrix<T>> t) {
		throw new ArithmeticException("Not a commutative ring!");
	}

	@Override
	public int krullDimension() {
		throw new ArithmeticException("Not a commutative ring!");
	}

	@Override
	public boolean isGeneratingAlgebra(List<Matrix<T>> s) {
		throw new ArithmeticException("Not a commutative ring!");
	}

	@Override
	public List<Matrix<T>> getAlgebraGenerators() {
		@SuppressWarnings("unchecked")
		T[][] matrix = (T[][]) Array.newInstance(zero().entry(1, 1).getClass(), dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				matrix[i][j] = ring.zero();
			}
		}
		List<Matrix<T>> generators = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			matrix[i][i] = ring.one();
			generators.add(new Matrix<T>(matrix));
			matrix[i][i] = ring.zero();
			if (i > 0) {
				matrix[i - 1][i] = ring.one();
				generators.add(new Matrix<T>(matrix));
				matrix[i - 1][i] = ring.zero();
				matrix[i][i - 1] = ring.one();
				generators.add(new Matrix<T>(matrix));
				matrix[i][i - 1] = ring.zero();
			}
		}
		return generators;
	}

	@Override
	public String toString() {
		return module.toString();
	}
}
