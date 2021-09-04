package fields.quaternions;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.quaternions.AbstractQuaternions.Quaternion;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.SubVectorSpace;
import fields.vectors.Vector;

public abstract class AbstractQuaternions<T extends Element<T>> extends AbstractAlgebra<T, Quaternion<T>>
		implements Quaternions<T> {
	public static class Quaternion<T extends Element<T>> extends AbstractElement<Quaternion<T>> {
		private T t;
		private T x;
		private T y;
		private T z;

		private Quaternion(T t, T x, T y, T z) {
			this.t = t;
			this.x = x;
			this.y = y;
			this.z = z;
		}

		public Vector<T> asVector() {
			return new Vector<>(t, x, y, z);
		}

		public T realPart() {
			return t;
		}

		public Vector<T> imaginaryPart() {
			return new Vector<>(x, y, z);
		}

		@Override
		public int compareTo(Quaternion<T> o) {
			int cmp = t.compareTo(o.t);
			if (cmp != 0) {
				return cmp;
			}
			cmp = x.compareTo(o.x);
			if (cmp != 0) {
				return cmp;
			}
			cmp = y.compareTo(o.y);
			if (cmp != 0) {
				return cmp;
			}
			return z.compareTo(o.z);
		}

		@Override
		public String toString() {
			return t + " + " + x + "*i + " + y + "*j + " + z + "*k";
		}

	}

	private final Field<T> base;
	private final T a;
	private final T b;
	private final FiniteVectorSpace<T> asVectorSpace;
	private List<Matrix<T>> matrixBasis;
	private final SubVectorSpace<T, Quaternion<T>> pureQuaternions;

	public AbstractQuaternions(Field<T> base, T a, T b) {
		if (a.equals(base.zero()) || b.equals(base.zero())) {
			throw new ArithmeticException("Ill defined");
		}
		this.base = base;
		this.a = a;
		this.b = b;
		this.asVectorSpace = new FiniteVectorSpace<>(base, 4);
		List<Quaternion<T>> pureQuaternions = new ArrayList<>();
		pureQuaternions.add(i());
		pureQuaternions.add(j());
		pureQuaternions.add(k());
		this.pureQuaternions = new SubVectorSpace<>(this, pureQuaternions);
	}

	@Override
	public String toString() {
		return "Q(" + base + ", i^2 = " + a + " j^2 = " + b + ")";
	}

	@Override
	public Quaternion<T> getEmbedding(T t) {
		return getElement(t, base.zero(), base.zero(), base.zero());
	}

	@Override
	public Quaternion<T> zero() {
		return getEmbedding(base.zero());
	}

	@Override
	public Quaternion<T> one() {
		return getEmbedding(base.one());
	}

	@Override
	public Quaternion<T> i() {
		return getElement(base.zero(), base.one(), base.zero(), base.zero());
	}

	@Override
	public Quaternion<T> j() {
		return getElement(base.zero(), base.zero(), base.one(), base.zero());
	}

	@Override
	public Quaternion<T> k() {
		return getElement(base.zero(), base.zero(), base.zero(), base.one());
	}

	@Override
	public Quaternion<T> getElement(T t, T x, T y, T z) {
		return new Quaternion<>(t, x, y, z);
	}

	@Override
	public T a() {
		return a;
	}

	@Override
	public T b() {
		return b;
	}

	@Override
	public T asRealPart(Quaternion<T> t) {
		if (!t.imaginaryPart().equals(new Vector<>(base.zero(), base.zero(), base.zero()))) {
			throw new ArithmeticException("Not an embedded field element");
		}
		return t.realPart();
	}

	@Override
	public SubVectorSpace<T, Quaternion<T>> pureQuaternions() {
		return pureQuaternions;
	}

	@Override
	public List<Quaternion<T>> getAlgebraGenerators() {
		List<Quaternion<T>> result = new ArrayList<>();
		result.add(i());
		result.add(j());
		return result;
	}

	@Override
	public boolean isGeneratingAlgebra(List<Quaternion<T>> s) {
		Quaternion<T> firstGenerator = null;
		for (Quaternion<T> t : s) {
			if (t.equals(conjugate(t))) {
				continue;
			}
			if (firstGenerator == null) {
				firstGenerator = t;
				continue;
			}
			if (!multiply(t, firstGenerator).equals(multiply(firstGenerator, t))) {
				return true;
			}
		}
		return false;
	}

	@Override
	public Quaternion<T> fromVector(Vector<T> t) {
		return getElement(t.get(1), t.get(2), t.get(3), t.get(4));
	}

	private List<Quaternion<T>> fromVectors(List<Vector<T>> ts) {
		List<Quaternion<T>> result = new ArrayList<>();
		for (Vector<T> t : ts) {
			result.add(fromVector(t));
		}
		return result;
	}

	@Override
	public Quaternion<T> add(Quaternion<T> s1, Quaternion<T> s2) {
		return getElement(base.add(s1.t, s2.t), base.add(s1.x, s2.x), base.add(s1.y, s2.y), base.add(s1.z, s2.z));
	}

	@Override
	public Quaternion<T> negative(Quaternion<T> s) {
		return getElement(base.negative(s.t), base.negative(s.x), base.negative(s.y), base.negative(s.z));
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public boolean isLinearIndependent(List<Quaternion<T>> s) {
		return asVectorSpace.isLinearIndependent(asVectors(s));
	}

	@Override
	public boolean isGeneratingModule(List<Quaternion<T>> s) {
		return asVectorSpace.isGeneratingModule(asVectors(s));
	}

	@Override
	public List<List<T>> nonTrivialCombinations(List<Quaternion<T>> s) {
		return asVectorSpace.nonTrivialCombinations(asVectors(s));
	}

	@Override
	public List<Quaternion<T>> getModuleGenerators() {
		return fromVectors(asVectorSpace.getModuleGenerators());
	}

	@Override
	public Vector<T> asVector(Quaternion<T> s) {
		return s.asVector();
	}

	private List<Vector<T>> asVectors(List<Quaternion<T>> ts) {
		List<Vector<T>> result = new ArrayList<>();
		for (Quaternion<T> t : ts) {
			result.add(asVector(t));
		}
		return result;
	}

	@Override
	public Exactness exactness() {
		return base.exactness();
	}

	@Override
	public Quaternion<T> getRandomElement() {
		return getElement(base.getRandomElement(), base.getRandomElement(), base.getRandomElement(),
				base.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return asVectorSpace.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return asVectorSpace.getNumberOfElements();
	}

	@Override
	public Iterator<Quaternion<T>> iterator() {
		return new Iterator<>() {
			private Iterator<Vector<T>> it = asVectorSpace.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public Quaternion<T> next() {
				return fromVector(it.next());
			}
		};
	}

	@Override
	public BigInteger characteristic() {
		return base.characteristic();
	}

	@Override
	public Quaternion<T> multiply(Quaternion<T> t1, Quaternion<T> t2) {
		return getElement(
				base.add(base.multiply(t1.t, t2.t), base.multiply(a, t1.x, t2.x), base.multiply(b, t1.y, t2.y),
						base.multiply(base.multiply(-1, a, b), t1.z, t2.z)),
				base.add(base.multiply(t1.t, t2.x), base.multiply(t1.x, t2.t), base.multiply(b, t1.z, t2.y),
						base.multiply(-1, b, t1.y, t2.z)),
				base.add(base.multiply(t1.t, t2.y), base.multiply(t1.y, t2.t), base.multiply(a, t1.x, t2.z),
						base.multiply(-1, a, t1.z, t2.x)),
				base.add(base.multiply(t1.t, t2.z), base.multiply(t1.z, t2.t), base.multiply(t1.x, t2.y),
						base.multiply(-1, t1.y, t2.x)));
	}

	@Override
	public boolean isUnit(Quaternion<T> t) {
		return !reducedNorm(t).equals(base.zero());
	}

	@Override
	public Quaternion<T> inverse(Quaternion<T> t) {
		return scalarMultiply(base.inverse(reducedNorm(t)), conjugate(t));
	}

	@Override
	public Quaternion<T> conjugate(Quaternion<T> t) {
		return getElement(t.t, base.negative(t.x), base.negative(t.y), base.negative(t.z));
	}

	@Override
	public T reducedNorm(Quaternion<T> t) {
		return base.add(base.multiply(t.t, t.t), base.multiply(-1, a, t.x, t.x), base.multiply(-1, b, t.y, t.y),
				base.multiply(base.multiply(a, b), t.z, t.z));
	}

	@Override
	public T reducedTrace(Quaternion<T> t) {
		return base.multiply(2, t.t);
	}

	@Override
	public QuadraticForm<T, Quaternion<T>> normForm() {
		PolynomialRing<T> polynomials = AbstractPolynomialRing.getPolynomialRing(base, 4, Monomial.GREVLEX);
		Map<Monomial, T> coeffs = new TreeMap<>();
		coeffs.put(polynomials.getMonomial(new int[] { 2, 0, 0, 0 }), base.one());
		coeffs.put(polynomials.getMonomial(new int[] { 0, 2, 0, 0 }), base.negative(a));
		coeffs.put(polynomials.getMonomial(new int[] { 0, 0, 2, 0 }), base.negative(b));
		coeffs.put(polynomials.getMonomial(new int[] { 0, 0, 0, 2 }), base.multiply(a, b));
		return new QuadraticForm<>(this, polynomials.getPolynomial(coeffs));
	}

	@Override
	public QuadraticForm<T, Quaternion<T>> restrictedNormForm() {
		PolynomialRing<T> polynomials = AbstractPolynomialRing.getPolynomialRing(base, 3, Monomial.GREVLEX);
		Map<Monomial, T> coeffs = new TreeMap<>();
		coeffs.put(polynomials.getMonomial(new int[] { 2, 0, 0 }), base.negative(a));
		coeffs.put(polynomials.getMonomial(new int[] { 0, 2, 0 }), base.negative(b));
		coeffs.put(polynomials.getMonomial(new int[] { 0, 0, 2 }), base.multiply(a, b));
		return new QuadraticForm<>(pureQuaternions(), polynomials.getPolynomial(coeffs));
	}

	@Override
	public boolean isZeroDivisor(Quaternion<T> t) {
		return reducedNorm(t).equals(base.zero());
	}

	@Override
	public abstract T discriminant();

	@Override
	public abstract int hilbertSymbol();

	@Override
	public abstract boolean isIntegral();

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isSplit() {
		return !isIntegral();
	}

	/**
	 * Returns a pure quaternion which has a square that is a square in the base
	 * field
	 * 
	 * @return
	 */
	protected abstract Quaternion<T> splittingElement();

	protected abstract Quaternion<T> normalize(Quaternion<T> t);

	protected Quaternion<T> negativeCommutor(Quaternion<T> t) {
		if (!t.realPart().equals(base.zero())) {
			throw new ArithmeticException("Not a pure Quaternion!");
		}
		if (t.equals(zero())) {
			throw new ArithmeticException("Zero!");
		}
		if (t.y.equals(base.zero()) && t.z.equals(base.zero())) {
			return j();
		}
		if (t.x.equals(base.zero()) && t.z.equals(base.zero())) {
			return k();
		}
		if (t.x.equals(base.zero()) && t.y.equals(base.zero())) {
			return i();
		}
		return normalize(getElement(base.zero(), base.multiply(b, t.z), base.multiply(a, t.z), base.add(t.x, t.y)));
	}

	private List<Matrix<T>> matrixBasisForSquare(T sqrtA, T b) {
		List<Matrix<T>> matrixBasis = new ArrayList<>();
		MatrixAlgebra<T> m = new FiniteVectorSpace<>(base, 2).matrixAlgebra();
		matrixBasis.add(m.one());
		List<Vector<T>> m1 = new ArrayList<>();
		m1.add(new Vector<>(sqrtA, base.zero()));
		m1.add(new Vector<>(base.zero(), base.negative(sqrtA)));
		Matrix<T> matrix1 = Matrix.fromRows(m1);
		matrixBasis.add(matrix1);
		List<Vector<T>> m2 = new ArrayList<>();
		m2.add(new Vector<>(base.zero(), b));
		m2.add(new Vector<>(base.one(), base.zero()));
		Matrix<T> matrix2 = Matrix.fromRows(m2);
		matrixBasis.add(matrix2);
		matrixBasis.add(m.multiply(matrix1, matrix2));
		return matrixBasis;
	}

	@Override
	public List<Matrix<T>> matrixBasis() {
		if (!isSplit()) {
			throw new ArithmeticException("Not a split quaternion algebra!");
		}
		if (matrixBasis == null) {
			MatrixAlgebra<T> m = new FiniteVectorSpace<>(base, 2).matrixAlgebra();
			if (base.hasSqrt(a())) {
				matrixBasis = matrixBasisForSquare(base.sqrt(a()).keySet().iterator().next(), b);
			} else {
				Quaternion<T> splittingElement = splittingElement();
				Quaternion<T> secondElement = negativeCommutor(splittingElement);
				T sqrtSplittingElement = base.sqrt(asRealPart(multiply(splittingElement, splittingElement))).keySet()
						.iterator().next();
				T secondElementSquare = asRealPart(multiply(secondElement, secondElement));
				List<Quaternion<T>> basis = new ArrayList<>();
				basis.add(one());
				basis.add(splittingElement);
				basis.add(secondElement);
				basis.add(multiply(splittingElement, secondElement));
				SubVectorSpace<T, Quaternion<T>> spaceInBasis = new SubVectorSpace<>(this, basis);
				List<Matrix<T>> matrixBasisForBasis = matrixBasisForSquare(sqrtSplittingElement, secondElementSquare);
				matrixBasis = new ArrayList<>();
				for (Quaternion<T> element : getBasis()) {
					Vector<T> asVector = spaceInBasis.asVector(element);
					Matrix<T> elementAsMatrix = m.zero();
					for (int i = 0; i < asVector.dimension(); i++) {
						elementAsMatrix = m.add(elementAsMatrix,
								m.scalarMultiply(asVector.get(i + 1), matrixBasisForBasis.get(i)));
					}
					matrixBasis.add(elementAsMatrix);
				}
			}
		}
		return matrixBasis;
	}

	@Override
	public Matrix<T> asMatrix(Quaternion<T> t) {
		MatrixAlgebra<T> m = new FiniteVectorSpace<>(base, 2).matrixAlgebra();
		Matrix<T> elementAsMatrix = m.zero();
		Vector<T> elementAsVector = asVector(t);
		for (int i = 0; i < 4; i++) {
			elementAsMatrix = m.add(elementAsMatrix,
					m.scalarMultiply(elementAsVector.get(i + 1), matrixBasis().get(i)));
		}
		return elementAsMatrix;
	}

	@Override
	public MathMap<Quaternion<T>, Matrix<T>> asMatrixMap() {
		return new MathMap<>() {

			@Override
			public Matrix<T> evaluate(Quaternion<T> t) {
				return asMatrix(t);
			}
		};
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public FactorizationResult<Quaternion<T>> uniqueFactorization(Quaternion<T> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isCommutative() {
		return false;
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isDedekindDomain() {
		return false;
	}

	@Override
	public boolean isDivisible(Quaternion<T> dividend, Quaternion<T> divisor) {
		return isUnit(divisor);
	}

	@Override
	public QuotientAndRemainderResult<Quaternion<T>> quotientAndRemainder(Quaternion<T> dividend,
			Quaternion<T> divisor) {
		if (!divisor.imaginaryPart().equals(new Vector<>(base.zero(), base.zero(), base.zero()))) {
			throw new ArithmeticException("Not commutative, use inverse and multiply instead!");
		}
		return new QuotientAndRemainderResult<>(scalarMultiply(base.inverse(divisor.realPart()), dividend), zero());
	}

	@Override
	public BigInteger euclidMeasure(Quaternion<T> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Quaternion<T> projectToUnit(Quaternion<T> t) {
		if (isZeroDivisor(t)) {
			return one();
		}
		return t;
	}

	@Override
	public Iterable<Quaternion<T>> getUnits() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int krullDimension() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public IdealResult<Quaternion<T>, ?> getIdealWithTransforms(List<Quaternion<T>> generators) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<Quaternion<T>> intersect(Ideal<Quaternion<T>> t1, Ideal<Quaternion<T>> t2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<Quaternion<T>> radical(Ideal<Quaternion<T>> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Field<T> getField() {
		return base;
	}

	@Override
	public Field<T> getRing() {
		return base;
	}

	@Override
	public List<Quaternion<T>> getBasis() {
		return getModuleGenerators();
	}

	@Override
	public int dimension() {
		return 4;
	}

	@Override
	public MatrixAlgebra<T> matrixAlgebra() {
		return asVectorSpace.matrixAlgebra();
	}

}
