package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractRing;
import fields.interfaces.Ideal;
import fields.interfaces.Ring;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;

public class ModuloNumberFieldIdeal extends AbstractRing<NFE> implements Ring<NFE> {
	private NumberFieldIntegers order;
	private NumberFieldIdeal ideal;

	public ModuloNumberFieldIdeal(NumberFieldIntegers order, NumberFieldIdeal ideal) {
		this.order = order;
		this.ideal = ideal;
	}

	public NFE getEmbedding(NFE t) {
		return ideal.residue(t);
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public NFE getRandomElement() {
		return getEmbedding(order.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return ideal.norm().getValue();
	}

	@Override
	public Iterator<NFE> iterator() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public NFE zero() {
		return order.zero();
	}

	@Override
	public NFE one() {
		return order.one();
	}

	@Override
	public BigInteger characteristic() {
		return ideal.intGenerator().getValue();
	}

	@Override
	public NFE add(NFE t1, NFE t2) {
		return getEmbedding(order.add(t1, t2));
	}

	@Override
	public NFE negative(NFE t) {
		return getEmbedding(order.negative(t));
	}

	@Override
	public NFE multiply(NFE t1, NFE t2) {
		return getEmbedding(order.multiply(t1, t2));
	}

	@Override
	public boolean isUnit(NFE t) {
		NumberFieldIdeal asIdeal = order.getIdeal(Collections.singletonList(t));
		return order.coprime(asIdeal, ideal);
	}

	@Override
	public NFE inverse(NFE t) {
		NumberFieldIdeal asIdeal = order.getIdeal(Collections.singletonList(t));
		BezoutIdentityResult<NFE> bezout = order.bezoutIdentity(asIdeal, ideal);
		return getEmbedding(bezout.getCoeff1().get(0));
	}

	@Override
	public boolean isCommutative() {
		return true;
	}

	@Override
	public boolean isIntegral() {
		return ideal.isMaximal();
	}

	@Override
	public boolean isReduced() {
		return order.idealFactorization(ideal).squareFree();
	}

	@Override
	public boolean isIrreducible() {
		return order.idealFactorization(ideal).primeFactors().size() == 1;
	}

	@Override
	public boolean isZeroDivisor(NFE t) {
		return !isUnit(t);
	}

	@Override
	public boolean isEuclidean() {
		return isIntegral();
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return isIntegral();
	}

	@Override
	public FactorizationResult<NFE, NFE> uniqueFactorization(NFE t) {
		if (!isIntegral()) {
			throw new ArithmeticException("Not a UFD");
		}
		return new FactorizationResult<>(t, Collections.emptySortedMap());
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return isIntegral();
	}

	@Override
	public boolean isDedekindDomain() {
		return isIntegral();
	}

	@Override
	public boolean isDivisible(NFE dividend, NFE divisor) {
		return order.add(ideal, order.getIdeal(Collections.singletonList(divisor))).contains(dividend);
	}

	@Override
	public QuotientAndRemainderResult<NFE> quotientAndRemainder(NFE dividend, NFE divisor) {
		List<NFE> generators = new ArrayList<>();
		generators.add(divisor);
		generators.addAll(ideal.generators());
		IdealResult<NFE, NumberFieldIdeal> asIdeal = order.getIdealWithTransforms(generators);
		List<NFE> generate = asIdeal.getIdeal().generate(dividend);
		NFE residue = getEmbedding(asIdeal.getIdeal().residue(dividend));
		NFE quotient = zero();
		for (int i = 0; i < generate.size(); i++) {
			quotient = add(quotient,
					getEmbedding(order.multiply(generate.get(i), asIdeal.getGeneratorExpressions().get(i).get(0))));
		}
		return new QuotientAndRemainderResult<>(quotient, residue);
	}

	@Override
	public BigInteger euclidMeasure(NFE t) {
		if (t.equals(zero())) {
		return null;
		}
		return BigInteger.ZERO;
	}

	@Override
	public NFE projectToUnit(NFE t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterable<NFE> getUnits() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int krullDimension() {
		return 0;
	}

	@Override
	public IdealResult<NFE, ?> getIdealWithTransforms(List<NFE> generators) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<NFE> intersect(Ideal<NFE> t1, Ideal<NFE> t2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<NFE> radical(Ideal<NFE> t) {
		// TODO Auto-generated method stub
		return null;
	}
}
