package fields.polynomials;

import java.math.BigInteger;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractField;
import fields.helper.FieldEmbedding;
import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.Element;
import fields.interfaces.FieldExtension;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing.ExtendedResultantResult;
import fields.local.LocalRingImplementation;
import fields.local.Value;
import util.Identity;

public class LocalizedUnivariatePolynomialRing<B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>>
		extends AbstractField<TExt<E>> implements DiscreteValuationField<TExt<E>, E> {
	private UnivariatePolynomial<E> uniformizer;
	private FE base;
	private TranscendentalFieldExtension<E> transcendental;
	private LocalRingImplementation<TExt<E>, E> integers;
	private FieldEmbedding<B, E, FE> reduction;
	private Reals r;

	public LocalizedUnivariatePolynomialRing(FE base, String variableName, UnivariatePolynomial<E> uniformizer) {
		this.base = base;
		this.uniformizer = uniformizer;
		this.transcendental = new TranscendentalFieldExtension<>(base, new String[] { variableName });
		this.integers = new LocalRingImplementation<>(this,
				base.getUnivariatePolynomialRing().toString() + "_(" + uniformizer + ")");
		this.reduction = base.getEmbeddedExtension(uniformizer);
		this.r = Reals.r(128);
	}

	@Override
	public Real value(TExt<E> t) {
		if (t.equals(transcendental.zero())) {
			return r.zero();
		}
		return r.exp(r.getInteger(-valuation(t).value()));
	}

	@Override
	public Reals getReals() {
		return r;
	}

	public TExt<E> getEmbedding(Polynomial<E> t) {
		return transcendental.getEmbedding(t);
	}

	@Override
	public TExt<E> zero() {
		return transcendental.zero();
	}

	@Override
	public TExt<E> one() {
		return transcendental.one();
	}

	@Override
	public BigInteger characteristic() {
		return transcendental.characteristic();
	}

	@Override
	public TExt<E> add(TExt<E> t1, TExt<E> t2) {
		return transcendental.add(t1, t2);
	}

	@Override
	public TExt<E> negative(TExt<E> t) {
		return transcendental.negative(t);
	}

	@Override
	public TExt<E> multiply(TExt<E> t1, TExt<E> t2) {
		return transcendental.multiply(t1, t2);
	}

	@Override
	public TExt<E> inverse(TExt<E> t) {
		return transcendental.inverse(t);
	}

	@Override
	public Exactness exactness() {
		return transcendental.exactness();
	}

	@Override
	public TExt<E> getRandomElement() {
		return transcendental.getRandomElement();
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
	public Iterator<TExt<E>> iterator() {
		return transcendental.iterator();
	}

	@Override
	public LocalizedUnivariatePolynomialRing<B, E, FE> asDedekindRing() {
		return this;
	}

	@Override
	public boolean isComplete() {
		return false;
	}

	@Override
	public TExt<E> inverse(TExt<E> t, int accuracy) {
		return inverse(t);
	}

	@Override
	public TExt<E> negative(TExt<E> t, int accuracy) {
		return negative(t);
	}

	@Override
	public Value valuation(TExt<E> t) {
		if (t.equals(transcendental.zero())) {
			return Value.INFINITY;
		}
		int value = 0;
		Polynomial<E> numerator = t.getNumerator();
		QuotientAndRemainderResult<Polynomial<E>> qr = base.getUnivariatePolynomialRing()
				.quotientAndRemainder(numerator, uniformizer);
		while (qr.getRemainder().equals(base.getUnivariatePolynomialRing().zero())) {
			numerator = qr.getQuotient();
			qr = base.getUnivariatePolynomialRing().quotientAndRemainder(numerator, uniformizer);
			value++;
		}
		Polynomial<E> denominator = t.getDenominator();
		qr = base.getUnivariatePolynomialRing().quotientAndRemainder(denominator, uniformizer);
		while (qr.getRemainder().equals(base.getUnivariatePolynomialRing().zero())) {
			denominator = qr.getQuotient();
			qr = base.getUnivariatePolynomialRing().quotientAndRemainder(denominator, uniformizer);
			value--;
		}
		return new Value(value);
	}

	@Override
	public TExt<E> uniformizer() {
		return transcendental.getEmbedding(uniformizer);
	}

	@Override
	public FE residueField() {
		return reduction.getField();
	}

	@Override
	public E reduceInteger(TExt<E> t) {
		return reduction.getField().divide(
				reduction.fromPolynomial(base.getUnivariatePolynomialRing().toUnivariate(t.getNumerator())),
				reduction.fromPolynomial(base.getUnivariatePolynomialRing().toUnivariate(t.getDenominator())));
	}

	@Override
	public TExt<E> upToUniformizerPower(TExt<E> t) {
		if (t.equals(transcendental.zero())) {
			transcendental.zero();
		}
		return transcendental.divide(t, transcendental.power(uniformizer(), valuation(t).value()));
	}

	@Override
	public TExt<E> liftToInteger(E s) {
		return transcendental.getEmbedding(reduction.asPolynomial(s));
	}

	@Override
	public TExt<E> round(TExt<E> t, int accuracy) {
		Value v = valuation(t);
		if (v.isInfinite()) {
			return t;
		}
		int add = 0;
		if (v.compareTo(Value.ZERO) < 0) {
			add = -v.value();
			t = multiply(t, power(uniformizer(), add));
		}
		List<Polynomial<E>> coefficients = base.getUnivariatePolynomialRing().adicDevelopment(t.getNumerator(),
				uniformizer);
		coefficients = coefficients.subList(0, Math.min(accuracy + add, coefficients.size()));
		Collections.reverse(coefficients);
		Polynomial<E> numerator = base.getUnivariatePolynomialRing().zero();
		for (Polynomial<E> c : coefficients) {
			numerator = base.getUnivariatePolynomialRing().multiply(uniformizer, numerator);
			numerator = base.getUnivariatePolynomialRing().add(numerator, c);
		}
		Polynomial<E> power = power(uniformizer(), accuracy + add).asInteger();
		ExtendedResultantResult<E> res = base.getUnivariatePolynomialRing().extendedResultant(t.getDenominator(),
				power);
		Polynomial<E> invertedDenominator = base.getUnivariatePolynomialRing().divideChecked(res.getCoeff1(),
				res.getGcd());
		numerator = base.getUnivariatePolynomialRing().multiply(numerator, invertedDenominator);
		return transcendental.getElement(numerator, base.getUnivariatePolynomialRing().power(uniformizer, add));
	}

	@Override
	public int getAccuracy() {
		throw new InfinityException();
	}

	@Override
	public LocalizedUnivariatePolynomialRing<B, E, FE> withAccuracy(int accuracy) {
		return this;
	}

	@Override
	public ExtensionOfDiscreteValuationField<TExt<E>, E, ?, ?, ?, ?, ?, ?, ?> getUniqueExtension(
			UnivariatePolynomial<TExt<E>> minimalPolynomial) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public OtherVersion<TExt<E>, TExt<E>, E, LocalizedUnivariatePolynomialRing<B, E, FE>> exact() {
		return new OtherVersion<>(this, new Identity<>(), new Identity<>());
	}

	@Override
	public OtherVersion<TExt<E>, ?, E, ?> complete(int accuracy) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TExt<E> getRandomInteger() {
		TExt<E> random = transcendental.getRandomElement();
		if (random.equals(transcendental.zero())) {
			return random;
		}
		Value v = valuation(random);
		return transcendental.multiply(
				transcendental.power(transcendental.getEmbedding(uniformizer), Math.max(0, -v.value())), random);
	}

	@Override
	public DiscreteValuationRing<TExt<E>, E> ringOfIntegers() {
		return integers;
	}
}
