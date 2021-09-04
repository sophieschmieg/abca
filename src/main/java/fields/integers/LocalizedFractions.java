package fields.integers;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractField;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.LocalField;
import fields.interfaces.LocalRing;
import fields.local.LocalRingImplementation;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.ValuationPivotStrategy;
import fields.local.Value;
import util.Identity;

public class LocalizedFractions extends AbstractField<Fraction> implements LocalField<Fraction, PFE> {
	private Fraction uniformizer;
	private Rationals q = Rationals.q();
	private Reals r = Reals.r(1024);
	private Integers z = Integers.z();
	private BigInteger prime;
	private PrimeField reduction;
	private LocalRing<Fraction, PFE> localRing;

	LocalizedFractions(BigInteger prime) {
		this.uniformizer = q.getEmbedding(prime);
		this.prime = prime;
		this.reduction = PrimeField.getPrimeField(prime);
		this.localRing = new LocalRingImplementation<>(this, "Z_(" + prime + ")");
	}

	@Override
	public boolean isComplete() {
		return false;
	}

	@Override
	public Real value(Fraction t) {
		if (t.equals(zero())) {
			return r.zero();
		}
		return r.exp(r.divide(r.log(r.getInteger(prime)), r.getInteger(valuation(t).value())));
	}

	@Override
	public Fraction zero() {
		return q.zero();
	}

	@Override
	public Fraction one() {
		return q.one();
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public Fraction add(Fraction t1, Fraction t2) {
		return q.add(t1, t2);
	}

	@Override
	public Fraction negative(Fraction t) {
		return q.negative(t);
	}

	@Override
	public Fraction multiply(Fraction t1, Fraction t2) {
		return q.multiply(t1, t2);
	}

	@Override
	public Fraction inverse(Fraction t) {
		return q.inverse(t);
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public Fraction getRandomElement() {
		return q.getRandomElement();
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
	public Iterator<Fraction> iterator() {
		return q.iterator();
	}

	@Override
	public Fraction inverse(Fraction t, int accuracy) {
		return inverse(t);
	}

	@Override
	public PivotStrategy<Fraction> preferredPivotStrategy() {
		return new ValuationPivotStrategy<>(this.valuation());
	}

	@Override
	public Fraction negative(Fraction t, int accuracy) {
		return negative(t);
	}

	@Override
	public boolean isInteger(Fraction t) {
		return valuation(t).compareTo(Value.ZERO) >= 0;
	}

	@Override
	public Value valuation(Fraction t) {
		if (t.equals(zero())) {
			return Value.INFINITY;
		}
		return new Value(z.valuation(t.getNumerator(), prime).value() - z.valuation(t.getDenominator(), prime).value());
	}

	@Override
	public Fraction uniformizer() {
		return uniformizer;
	}

	@Override
	public PrimeField reduction() {
		return reduction;
	}

	@Override
	public PFE reduce(Fraction t) {
		t.canonicalize();
		return reduction.divide(z.reduce(t.getNumerator(), prime), z.reduce(t.getDenominator(), prime));
	}

	@Override
	public Fraction upToUniformizerPower(Fraction t) {
		if (t.equals(zero())) {
			return zero();
		}
		return divide(t, power(uniformizer(), valuation(t).value()));
	}

	@Override
	public Fraction lift(PFE s) {
		return q.getEmbedding(z.lift(s));
	}

	@Override
	public Fraction round(Fraction t, int accuracy) {
		int denominatorValue = 0;
		BigInteger denominator = t.getDenominator().getValue();
		while (denominator.mod(prime).equals(BigInteger.ZERO)) {
			denominator = denominator.divide(prime);
			denominatorValue++;
		}
		if (denominatorValue + accuracy == 0) {
			return zero();
		}
		BigInteger mod = prime.pow(denominatorValue + accuracy);
		BigInteger denominatorInverted = denominator.modInverse(mod);
		Fraction numerator = getInteger(t.getNumerator().getValue().multiply(denominatorInverted).mod(mod));
		Fraction denominatorFraction = getInteger(prime.pow(denominatorValue));
		return divide(numerator, denominatorFraction);
	}

	@Override
	public int getAccuracy() {
		throw new InfinityException();
	}

	@Override
	public LocalField<Fraction, PFE> withAccuracy(int accuracy) {
		return this;
	}

	@Override
	public OtherVersion<Fraction, Fraction, PFE, LocalizedFractions> exact() {
		return new OtherVersion<>(this, new Identity<>(), new Identity<>());
	}

	@Override
	public OtherVersion<Fraction, PAdicNumber, PFE, PAdicField> complete(int accuracy) {
		PAdicField field = new PAdicField(prime, accuracy);
		return new OtherVersion<>(field, field.fromRationalMap(), field.toRationalMap());
	}

	@Override
	public Fraction getRandomInteger() {
		IntE numerator = z.getRandomElement();
		IntE denominator;
		do {
			denominator = z.getRandomElement();
		} while (denominator.equals(z.zero()));

		IntE primeInt = z.getInteger(prime);
		while (z.isDivisible(denominator, primeInt)) {
			denominator = z.divideChecked(denominator, primeInt);
		}
		return q.getFraction(numerator, denominator);
	}

	@Override
	public LocalRing<Fraction, PFE> ringOfIntegers() {
		return localRing;
	}

}
