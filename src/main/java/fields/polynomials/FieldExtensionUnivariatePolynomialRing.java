package fields.polynomials;

import java.io.IOException;
import java.math.BigInteger;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractRing;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DedekindRing;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.Element;
import fields.interfaces.FieldExtension;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.Value;
import util.FunctionMathMap;
import util.PeekableReader;

public class FieldExtensionUnivariatePolynomialRing<B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>>
		extends AbstractRing<Polynomial<E>> implements DedekindRing<Polynomial<E>, TExt<E>, E> {
	private FE base;
	private UnivariateTranscendentalFieldExtension<B, E, FE> quotientField;
	private UnivariatePolynomialRing<E> asUnivariatePolynomialRing;

	public FieldExtensionUnivariatePolynomialRing(FE base) {
		this(base, "X");
	}

	public FieldExtensionUnivariatePolynomialRing(FE base, String variableName) {
		super(true);
		this.base = base;
		this.asUnivariatePolynomialRing = base.getUnivariatePolynomialRing();
	}

	public FE getField() {
		return base;
	}

	@Override
	public FieldExtensionUnivariatePolynomialRing<B, E, FE> asDedekindRing() {
		return this;
	}

	@Override
	public FieldOfFractionsResult<Polynomial<E>, TExt<E>> fieldOfFractions() {
		return new FieldOfFractionsResult<>(this, quotientField(),
				new FunctionMathMap<>((Polynomial<E> t) -> quotientField().getEmbedding(t)),
				new FunctionMathMap<>((TExt<E> t) -> t.getNumerator()),
				new FunctionMathMap<>((TExt<E> t) -> t.getDenominator()),
				new FunctionMathMap<>((TExt<E> t) -> t.asInteger()));
	}

	@Override
	public Value valuation(Polynomial<E> t, Ideal<Polynomial<E>> maximalIdeal) {
		if (t.equals(zero())) {
			return Value.INFINITY;
		}
		Polynomial<E> generator = maximalIdeal.generators().get(0);
		int order = 0;
		QuotientAndRemainderResult<Polynomial<E>> qr = quotientAndRemainder(t, generator);
		while (qr.getRemainder().equals(zero())) {
			order++;
			qr = quotientAndRemainder(qr.getQuotient(), generator);
		}
		return new Value(order);
	}

	@Override
	public boolean isInteger(TExt<E> t) {
		return isUnit(t.getDenominator());
	}

	@Override
	public Polynomial<E> asInteger(TExt<E> t) {
		return t.asInteger();
	}

	@Override
	public Polynomial<E> getDenominator(TExt<E> t) {
		return t.getDenominator();
	}

	@Override
	public FE reduction(Ideal<Polynomial<E>> maximalIdeal) {
		return base.getEmbeddedExtension(asUnivariatePolynomialRing.toUnivariate(maximalIdeal.generators().get(0)))
				.getField();
	}

	@Override
	public E reduce(Polynomial<E> t, Ideal<Polynomial<E>> maximalIdeal) {
		return base.getEmbeddedExtension(asUnivariatePolynomialRing.toUnivariate(maximalIdeal.generators().get(0)))
				.fromPolynomial(asUnivariatePolynomialRing.toUnivariate(t));
	}

	@Override
	public Polynomial<E> lift(E s, Ideal<Polynomial<E>> maximalIdeal) {
		return base.getEmbeddedExtension(asUnivariatePolynomialRing.toUnivariate(maximalIdeal.generators().get(0)))
				.asPolynomial(s);
	}

	@Override
	public DiscreteValuationRing<TExt<E>, E> localize(Ideal<Polynomial<E>> maximalIdeal) {
		return localizeAndQuotient(maximalIdeal).ringOfIntegers();
	}

	@Override
	public DiscreteValuationField<TExt<E>, E> localizeAndQuotient(Ideal<Polynomial<E>> maximalIdeal) {
		return new LocalizedUnivariatePolynomialRing<>(base, asUnivariatePolynomialRing.getVariableName(),
				base.getUnivariatePolynomialRing().toUnivariate(maximalIdeal.principalGenerator()));
	}

	@Override
	public UnivariateTranscendentalFieldExtension<B, E, FE> quotientField() {
		if (quotientField == null) {
			quotientField = new UnivariateTranscendentalFieldExtension<>(base);
		}
		return quotientField;
	}

	@Override
	public Exactness exactness() {
		return base.exactness();
	}

	@Override
	public UnivariatePolynomial<E> parse(PeekableReader reader) throws IOException {
		return asUnivariatePolynomialRing.parse(reader);
	}

	@Override
	public UnivariatePolynomial<E> zero() {
		return asUnivariatePolynomialRing.zero();
	}

	@Override
	public UnivariatePolynomial<E> one() {
		return asUnivariatePolynomialRing.one();
	}

	@Override
	public Polynomial<E> getRandomElement() {
		return asUnivariatePolynomialRing.getRandomElement();
	}

	@Override
	public UnivariatePolynomial<E> add(Polynomial<E> s1, Polynomial<E> s2) {
		return asUnivariatePolynomialRing.add(s1, s2);
	}

	@Override
	public UnivariatePolynomial<E> negative(Polynomial<E> s) {
		return asUnivariatePolynomialRing.negative(s);
	}

	@Override
	public UnivariatePolynomial<E> multiply(Polynomial<E> t1, Polynomial<E> t2) {
		return asUnivariatePolynomialRing.multiply(t1, t2);
	}

	@Override
	public boolean isZeroDivisor(Polynomial<E> t) {
		return false;
	}

	@Override
	public boolean isUnit(Polynomial<E> t) {
		return t.degree() == 0;
	}

	@Override
	public UnivariatePolynomial<E> inverse(Polynomial<E> t) {
		return asUnivariatePolynomialRing.inverse(t);
	}

	@Override
	public boolean isDivisible(Polynomial<E> dividend, Polynomial<E> divisor) {
		return asUnivariatePolynomialRing.isDivisible(dividend, divisor);
	}

	@Override
	public BigInteger euclidMeasure(Polynomial<E> t) {
		return BigInteger.valueOf(t.degree());
	}

	@Override
	public Iterable<Polynomial<E>> getUnits() {
		return asUnivariatePolynomialRing.getUnits();
	}

	@Override
	public Iterator<Polynomial<E>> iterator() {
		return asUnivariatePolynomialRing.iterator();
	}

	@Override
	public UnivariatePolynomial<E> gcd(Polynomial<E> t1, Polynomial<E> t2) {
		return asUnivariatePolynomialRing.gcd(t1, t2);
	}

	@Override
	public ExtendedEuclideanResult<Polynomial<E>> extendedEuclidean(Polynomial<E> t1, Polynomial<E> t2) {
		return asUnivariatePolynomialRing.extendedEuclidean(t1, t2);
	}

	@Override
	public BigInteger characteristic() {
		return asUnivariatePolynomialRing.characteristic();
	}

	@Override
	public boolean isCommutative() {
		return true;
	}

	@Override
	public boolean isEuclidean() {
		return true;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return true;
	}

	@Override
	public FactorizationResult<Polynomial<E>, Polynomial<E>> uniqueFactorization(Polynomial<E> t) {
		return asUnivariatePolynomialRing.uniqueFactorization(t);
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return true;
	}

	@Override
	public QuotientAndRemainderResult<Polynomial<E>> quotientAndRemainder(Polynomial<E> dividend,
			Polynomial<E> divisor) {
		return asUnivariatePolynomialRing.quotientAndRemainder(dividend, divisor);
	}

	@Override
	public Polynomial<E> projectToUnit(Polynomial<E> t) {
		return asUnivariatePolynomialRing.projectToUnit(t);
	}

	@Override
	public int krullDimension() {
		return 1;
	}

	@Override
	public IdealResult<Polynomial<E>, PolynomialIdeal<E>> getIdealWithTransforms(List<Polynomial<E>> generators) {
		return asUnivariatePolynomialRing.getIdealWithTransforms(generators);
	}

	@Override
	public PolynomialIdeal<E> intersect(Ideal<Polynomial<E>> t1, Ideal<Polynomial<E>> t2) {
		return asUnivariatePolynomialRing.intersect(t1, t2);
	}

	@Override
	public PolynomialIdeal<E> radical(Ideal<Polynomial<E>> t) {
		return asUnivariatePolynomialRing.radical(t);
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return asUnivariatePolynomialRing.getNumberOfElements();
	}

}
