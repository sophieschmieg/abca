package varieties.curves.elliptic;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberFieldIntegers;
import fields.quaternions.RationalQuaternionOrder;
import fields.vectors.Vector;
import varieties.curves.elliptic.EllipticCurve.Isomorphism;
import varieties.curves.elliptic.EndomorphismRing.Endomorphism;

public class EndomorphismRing<T extends Element<T>> extends AbstractAlgebra<IntE, Endomorphism<T>> {
	public static class Endomorphism<T extends Element<T>> extends AbstractElement<Endomorphism<T>> {
		private Vector<IntE> asVector;

		private Endomorphism(Vector<IntE> asVector) {
			this.asVector = asVector;
		}

		@Override
		public int compareTo(Endomorphism<T> o) {
			return asVector.compareTo(o.asVector);
		}
	}

	private EllipticCurve<T> curve;
	private List<Isogeny<T>> basis;
	private Integers z;
	private NumberFieldIntegers asNumberFieldOrder;
	private RationalQuaternionOrder asQuaternionOrder;

	EndomorphismRing(EllipticCurve<T> curve) {
		this.curve = curve;
		this.z = Integers.z();
		Isomorphism<T> isomorphism = curve.identity();
		basis.add(isomorphism);
		if (!curve.getField().isFinite()) {
			throw new UnsupportedOperationException("Only curves over finite fields currently supported");
		}
		Isogeny<T> frobenius = curve.frobeniousEndomorphism();
		basis.add(frobenius);
		if (curve.isSupersingular()) {

		} else {
			Rationals q = Rationals.q();
			Fraction trace = q.getInteger(curve.trace());
			Fraction norm = q.getInteger(frobenius.getDegree());
			UnivariatePolynomialRing<Fraction> rationalPolynomials = q.getUnivariatePolynomialRing();
			NumberField nf = new NumberField(rationalPolynomials.getPolynomial(norm, trace, q.one()));
			asNumberFieldOrder = nf.maximalOrder();
		}
	}

	@Override
	public Endomorphism<T> getEmbedding(IntE t) {
		List<IntE> list = new ArrayList<>();
		list.add(t);
		for (int i = 1; i < basis.size(); i++) {
			list.add(z.zero());
		}
		return new Endomorphism<>(new Vector<>(list));
	}

	public IntE reducedNorm(Isogeny<T> isogeny) {
		return z.getInteger(isogeny.getDegree());
	}

	public IntE reducedTrace(Isogeny<T> isogeny) {
		return z.getInteger(curve.asMultiplicationIsogenyModOrder(new AdditionIsogeny<>(isogeny, isogeny.getDual())));
	}

	@Override
	public boolean isGeneratingAlgebra(List<Endomorphism<T>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public List<Endomorphism<T>> getAlgebraGenerators() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ring<IntE> getRing() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Endomorphism<T> zero() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Endomorphism<T> add(Endomorphism<T> s1, Endomorphism<T> s2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Endomorphism<T> negative(Endomorphism<T> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isFree() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Ideal<IntE> annihilator() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isLinearIndependent(List<Endomorphism<T>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isGeneratingModule(List<Endomorphism<T>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public List<List<IntE>> nonTrivialCombinations(List<Endomorphism<T>> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Endomorphism<T>> getModuleGenerators() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Vector<IntE> asVector(Endomorphism<T> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Exactness exactness() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Endomorphism<T> getRandomElement() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isFinite() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterator<Endomorphism<T>> iterator() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Endomorphism<T> one() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public BigInteger characteristic() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Endomorphism<T> multiply(Endomorphism<T> t1, Endomorphism<T> t2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isUnit(Endomorphism<T> t) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Endomorphism<T> inverse(Endomorphism<T> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isCommutative() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isIntegral() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isReduced() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isIrreducible() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isZeroDivisor(Endomorphism<T> t) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isEuclidean() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public FactorizationResult<Endomorphism<T>, Endomorphism<T>> uniqueFactorization(Endomorphism<T> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isDedekindDomain() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isDivisible(Endomorphism<T> dividend, Endomorphism<T> divisor) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public QuotientAndRemainderResult<Endomorphism<T>> quotientAndRemainder(Endomorphism<T> dividend,
			Endomorphism<T> divisor) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public BigInteger euclidMeasure(Endomorphism<T> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Endomorphism<T> projectToUnit(Endomorphism<T> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterable<Endomorphism<T>> getUnits() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int krullDimension() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public IdealResult<Endomorphism<T>, ?> getIdealWithTransforms(List<Endomorphism<T>> generators) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<Endomorphism<T>> intersect(Ideal<Endomorphism<T>> t1, Ideal<Endomorphism<T>> t2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<Endomorphism<T>> radical(Ideal<Endomorphism<T>> t) {
		// TODO Auto-generated method stub
		return null;
	}
}
