package fields.quaternions;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.interfaces.DedekindRing;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.quaternions.AbstractQuaternions.Quaternion;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

public abstract class AbstractQuaternionOrder<T extends Element<T>, I extends Element<I>, R extends Element<R>>
		extends AbstractAlgebra<I, Quaternion<T>> implements QuaternionOrder<T, I, R> {
	private DedekindRing<I, T, R> ring;
	private FieldOfFractionsResult<I, T> fieldOfFractions;
	private Field<T> field;
	private Quaternions<T> quaternions;
	private FreeModule<I> asFreeModule;
	private MatrixAlgebra<T> matrixAlgebra;
	private Matrix<T> toIntegralBasisBaseChange;
	private Matrix<T> fromIntegralBasisBaseChange;
	private List<Quaternion<T>> basis;

	public AbstractQuaternionOrder(DedekindRing<I, T, R> ring, Quaternions<T> quaternions, Quaternion<T> t1,
			Quaternion<T> t2, Quaternion<T> t3) {
		this.ring = ring;
		this.quaternions = quaternions;
		this.fieldOfFractions = ring.fieldOfFractions();
		this.field = fieldOfFractions.getField();
		if (!this.field.equals(quaternions.getField())) {
			throw new ArithmeticException("Fields do not match!");
		}
		this.basis = new ArrayList<>();
		this.basis.add(one());
		this.basis.add(t1);
		this.basis.add(t2);
		this.basis.add(t3);
		for (Quaternion<T> basisVector : basis) {
			if (!ring.isInteger(quaternions.reducedTrace(basisVector))
					|| !ring.isInteger(quaternions.reducedNorm(basisVector))) {
				throw new ArithmeticException("Not an integral basis!");
			}
		}
		List<Vector<T>> basisAsVectors = new ArrayList<>();
		basisAsVectors.add(quaternions.asVector(one()));
		basisAsVectors.add(quaternions.asVector(t1));
		basisAsVectors.add(quaternions.asVector(t2));
		basisAsVectors.add(quaternions.asVector(t3));
		this.asFreeModule = new FreeModule<>(ring, 4);
		this.matrixAlgebra = quaternions.matrixAlgebra();
		this.fromIntegralBasisBaseChange = Matrix.fromColumns(basisAsVectors);
		this.toIntegralBasisBaseChange = matrixAlgebra.inverse(fromIntegralBasisBaseChange);
	}

	@Override
	public Quaternion<T> getEmbedding(I t) {
		return quaternions.getEmbedding(fieldOfFractions.getEmbedding().evaluate(t));
	}

	@Override
	public boolean isGeneratingAlgebra(List<Quaternion<T>> s) {
		List<Quaternion<T>> asModuleGenerators = new ArrayList<>();
		asModuleGenerators.add(one());
		asModuleGenerators.addAll(s);
		for (Quaternion<T> t1 : s) {
			for (Quaternion<T> t2 : s) {
				asModuleGenerators.add(multiply(t1, t2));
			}
		}
		return isGeneratingModule(asModuleGenerators);
	}

	@Override
	public List<Quaternion<T>> getAlgebraGenerators() {
		return basis.subList(1, 4);
	}

	@Override
	public DedekindRing<I, T, R> getRing() {
		return ring;
	}

	@Override
	public Quaternion<T> zero() {
		return quaternions.zero();
	}

	@Override
	public Quaternion<T> add(Quaternion<T> s1, Quaternion<T> s2) {
		return quaternions.add(s1, s2);
	}

	@Override
	public Quaternion<T> negative(Quaternion<T> s) {
		return quaternions.negative(s);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public boolean isLinearIndependent(List<Quaternion<T>> s) {
		return asFreeModule.isLinearIndependent(asVectors(s));
	}

	@Override
	public boolean isGeneratingModule(List<Quaternion<T>> s) {
		return asFreeModule.isGeneratingModule(asVectors(s));
	}

	@Override
	public List<List<I>> nonTrivialCombinations(List<Quaternion<T>> s) {
		return asFreeModule.nonTrivialCombinations(asVectors(s));
	}

	@Override
	public List<Quaternion<T>> getModuleGenerators() {
		return this.basis;
	}

	@Override
	public Vector<I> asVector(Quaternion<T> s) {
		Vector<T> asVector = matrixAlgebra.multiply(toIntegralBasisBaseChange, quaternions.asVector(s));
		List<I> asIntegers = new ArrayList<>();
		for (T t : asVector.asList()) {
			asIntegers.add(ring.asInteger(t));
		}
		return new Vector<>(asIntegers);
	}

	private List<Vector<I>> asVectors(List<Quaternion<T>> s) {
		List<Vector<I>> result = new ArrayList<>();
		for (Quaternion<T> t : s) {
			result.add(asVector(t));
		}
		return result;
	}

	@Override
	public boolean isElement(Quaternion<T> t) {
		Vector<T> asVector = matrixAlgebra.multiply(toIntegralBasisBaseChange, quaternions.asVector(t));
		for (T e : asVector.asList()) {
			if (!ring.isInteger(e)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public Quaternion<T> getRandomElement() {
		return fromVector(asFreeModule.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<Quaternion<T>> iterator() {
		throw new InfinityException();
	}

	@Override
	public Quaternion<T> one() {
		return quaternions.one();
	}

	@Override
	public BigInteger characteristic() {
		return quaternions.characteristic();
	}

	@Override
	public Quaternion<T> multiply(Quaternion<T> t1, Quaternion<T> t2) {
		return quaternions.multiply(t1, t2);
	}

	@Override
	public boolean isUnit(Quaternion<T> t) {
		return isElement(quaternions.inverse(t));
	}

	@Override
	public Quaternion<T> inverse(Quaternion<T> t) {
		Quaternion<T> inverse = quaternions.inverse(t);
		if (!isElement(inverse)) {
			throw new ArithmeticException("Not a unit");
		}
		return inverse;
	}

	@Override
	public boolean isCommutative() {
		return false;
	}

	@Override
	public boolean isIntegral() {
		return quaternions.isIntegral();
	}
	
	@Override
	public boolean isReduced() {
		return quaternions.isReduced();
	}
	
	@Override
	public boolean isIrreducible() {
		return quaternions.isIrreducible();
	}

	@Override
	public boolean isZeroDivisor(Quaternion<T> t) {
		return quaternions.isZeroDivisor(t);
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
	public FactorizationResult<Quaternion<T>, Quaternion<T>> uniqueFactorization(Quaternion<T> t) {
		throw new ArithmeticException("Not a UFD!");
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
	public boolean isDivisible(Quaternion<T> dividend, Quaternion<T> divisor) {
		Quaternion<T> quotient = quaternions.divide(dividend, divisor);
		return isElement(quotient);
	}

	@Override
	public Quaternion<T> divideChecked(Quaternion<T> dividend, Quaternion<T> divisor) {
		Quaternion<T> quotient = quaternions.divide(dividend, divisor);
		if (!isElement(quotient)) {
			throw new ArithmeticException("Not divisible!");
		}
		return quotient;
	}

	@Override
	public QuotientAndRemainderResult<Quaternion<T>> quotientAndRemainder(Quaternion<T> dividend,
			Quaternion<T> divisor) {
		Quaternion<T> quotient = quaternions.divide(dividend, divisor);
		if (!isElement(quotient)) {
			return new QuotientAndRemainderResult<>(zero(), dividend);
		}
		return new QuotientAndRemainderResult<>(quotient, zero());
	}

	@Override
	public BigInteger euclidMeasure(Quaternion<T> t) {
		throw new ArithmeticException("not an euclidean ring");
	}

	@Override
	public Quaternion<T> projectToUnit(Quaternion<T> t) {
		return one();
	}

	@Override
	public Iterable<Quaternion<T>> getUnits() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int krullDimension() {
		throw new ArithmeticException("not a commutative ring!");
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
	public Quaternions<T> getQuaternions() {
		return quaternions;
	}

	@Override
	public I reducedDiscriminant() {
		Quaternion<T> b12 = quaternions.multiply(basis.get(1), basis.get(2));
		Quaternion<T> b21 = quaternions.multiply(basis.get(2), basis.get(1));
		Quaternion<T> m = quaternions.multiply(quaternions.subtract(b12, b21), quaternions.conjugate(basis.get(3)));
		return ring.asInteger(quaternions.reducedTrace(m));
	}

	@Override
	public boolean isMaximal() {
		return quaternions.discriminant().equals(fieldOfFractions.getEmbedding().evaluate(reducedDiscriminant()));
	}
}
