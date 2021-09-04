package fields.local;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.Random;

import fields.exceptions.InfinityException;
import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.ModularIntegerRing;
import fields.finitefields.ModularIntegerRing.ModularIntegerRingElement;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.LocalizedFractions;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Field;
import fields.interfaces.LocalField;
import fields.interfaces.LocalRing;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.local.CompleteLocalFieldExtension.Ext;
import fields.local.PAdicField.PAdicNumber;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.ValuationPivotStrategy;
import util.Identity;

public class PAdicField extends AbstractField<PAdicNumber> implements Field<PAdicNumber>, LocalField<PAdicNumber, PFE> {
	private BigInteger prime;
	private PrimeField reduced;
	private int accuracy;
	private BigInteger modulus;
	private LocalRing<PAdicNumber, PFE> localRing;
	private Reals r = Reals.r(1024);

	public class PAdicNumber extends AbstractElement<PAdicNumber> {
		private BigInteger value;
		private int lowestPower;

		private PAdicNumber(int lowestPower, BigInteger value) {
			if (lowestPower > accuracy) {
				value = BigInteger.ZERO;
				lowestPower = 0;
			}
			while (!value.equals(BigInteger.ZERO) && value.mod(prime).equals(BigInteger.ZERO)) {
				value = value.divide(prime);
				lowestPower++;
			}
			BigInteger modulus = PAdicField.this.modulus;
			if (lowestPower < 0) {
				modulus.multiply(prime.pow(-lowestPower));
			} /*
				 * else if (lowestPower > 0) { modulus.divide(prime.pow(lowestPower)); }
				 */
			value = value.mod(modulus);
			if (value.equals(BigInteger.ZERO)) {
				lowestPower = 0;
			}
			if (lowestPower >= accuracy) {
				lowestPower = 0;
				value = BigInteger.ZERO;
			}
			this.lowestPower = lowestPower;
			this.value = value;
		}

		public String toString() {
			StringBuilder build = new StringBuilder().append("...");
			boolean first = true;
			for (int power = getAccuracy() - 1; power >= Math.min(lowestPower, 0); power--) {
				if (prime.compareTo(BigInteger.TEN) > 0 && !first) {
					build.append(" ");
				}
				if (power == -1) {
					build.append(".");
				}
				build.append(digit(power));
				first = false;
			}
			return build.toString();
		}

		public BigInteger digit(int position) {
			position -= lowestPower;
			if (position < 0) {
				return BigInteger.ZERO;
			}
			return value.divide(prime.pow(position)).mod(prime);
		}

		@Override
		public int compareTo(PAdicNumber o) {
			PAdicNumber roundedThis = round(this, getAccuracy());
			PAdicNumber roundedOther = round(o, getAccuracy());
			boolean thisZero = false;
			if (roundedThis.value.equals(BigInteger.ZERO)) {
				thisZero = true;
			}
			boolean otherZero = false;
			if (roundedOther.value.equals(BigInteger.ZERO)) {
				otherZero = true;
			}
			if (thisZero && otherZero) {
				return 0;
			}
			if (thisZero) {
				return 1;
			}
			if (otherZero) {
				return -1;
			}
			if (roundedThis.lowestPower != roundedOther.lowestPower) {
				return roundedThis.lowestPower - roundedOther.lowestPower;
			}
			return roundedThis.value.compareTo(roundedOther.value);
		}

	}

	public PAdicField(BigInteger prime, int accuracy) {
		if (!prime.isProbablePrime(100)) {
			throw new ArithmeticException("not a prime!");
		}
		accuracy = Math.max(accuracy, 3);
		this.accuracy = accuracy;
		this.prime = prime;
		this.modulus = this.prime.pow(accuracy);
		this.reduced = PrimeField.getPrimeField(prime);
		this.localRing = new LocalRingImplementation<>(this, "Z_" + prime);
	}

	@Override
	public Exactness exactness() {
		return Exactness.FIXED_POINT;
	}

	@Override
	public boolean isComplete() {
		return true;
	}

	@Override
	public Extension<PAdicNumber, PAdicNumber, Ext<PAdicNumber>, CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField>> getExtension(
			UnivariatePolynomial<PAdicNumber> minimalPolynomial) {
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> extension = new CompleteLocalFieldExtension<>(
				minimalPolynomial, this, FiniteField.getFiniteField(prime));
		return new Extension<>(extension, this, extension.getEmbeddingMap(), extension.asVectorMap());
	}

	@Override
	public int getAccuracy() {
		return accuracy;
	}

	@Override
	public PAdicField withAccuracy(int accuracy) {
		return new PAdicField(prime, accuracy);
	}

	@Override
	public OtherVersion<PAdicNumber, Fraction, PFE, LocalizedFractions> exact() {
		return new OtherVersion<>(Rationals.q().withValuation(prime), toRationalMap(), fromRationalMap());
	}

	@Override
	public OtherVersion<PAdicNumber, PAdicNumber, PFE, PAdicField> complete(int accuracy) {
		return new OtherVersion<>(withAccuracy(accuracy), new Identity<>(), new Identity<>());
	}

	public PrimeField reduction() {
		return reduced;
	}
	
	public ModularIntegerRing reduction(int numDigits) {
		return new ModularIntegerRing(prime.pow(numDigits));
	}

	@Override
	public PAdicNumber getRandomElement() {
		Random r = new Random();
		BigInteger powerInt = new BigInteger(modulus.bitLength(), r);
		PAdicNumber powerNumber = new PAdicNumber(0, powerInt);
		int lowestPower = -powerNumber.lowestPower;
		BigInteger value = new BigInteger(prime.pow(accuracy - lowestPower).bitLength(), r);
		return new PAdicNumber(lowestPower, value);
	}

	@Override
	public PAdicNumber getRandomInteger() {
		Random r = new Random();
		BigInteger value = new BigInteger(modulus.bitLength(), r);
		return new PAdicNumber(0, value);
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
	public Iterator<PAdicNumber> iterator() {
		throw new InfinityException();
	}

	@Override
	public PAdicNumber zero() {
		return new PAdicNumber(0, BigInteger.ZERO);
	}

	@Override
	public PAdicNumber one() {
		return new PAdicNumber(0, BigInteger.ONE);
	}

	@Override
	public Real value(PAdicNumber t) {
		if (t.equals(zero())) {
			return r.zero();
		}
		return r.exp(r.divide(r.log(r.getInteger(prime)), r.getInteger(valuation(t).value())));

	}

	@Override
	public Value valuation(PAdicNumber p) {
		if (p.value.equals(BigInteger.ZERO)) {
			return Value.INFINITY;
		}
		return new Value(p.lowestPower);
	}

	public Fraction toRational(PAdicNumber t) {
		Rationals q = Rationals.q();
		if (t.lowestPower >= accuracy) {
			return q.zero();
		}
		BigInteger modulus = prime.pow(accuracy - t.lowestPower - 1);
		BigInteger period = t.value.divide(modulus);
		BigInteger denominator = BigInteger.ONE.subtract(prime);
		BigInteger numerator = period.multiply(prime.pow(accuracy - 1));
		Integers z = Integers.z();
		return q.add(q.getFraction(z.getInteger(numerator), z.getInteger(denominator)),
				q.getFraction(z.getInteger(t.value.mod(modulus).multiply(prime.pow(Math.max(0, t.lowestPower)))),
						z.getInteger(prime.pow(Math.max(0, -t.lowestPower)))));
	}

	public MathMap<PAdicNumber, Fraction> toRationalMap() {
		return new MathMap<>() {
			@Override
			public Fraction evaluate(PAdicNumber n) {
				return toRational(n);
			}
		};
	}

	public PAdicNumber fromRational(Fraction t) {
		BigInteger numerator = t.getNumerator().getValue();
		BigInteger denominator = t.getDenominator().getValue();
		return divide(getInteger(numerator), getInteger(denominator));
	}

	public MathMap<Fraction, PAdicNumber> fromRationalMap() {
		return new MathMap<>() {
			@Override
			public PAdicNumber evaluate(Fraction n) {
				return fromRational(n);
			}
		};
	}

	@Override
	public PFE reduce(PAdicNumber p) {
		PAdicNumber t = (PAdicNumber) p;
		if (t.lowestPower < 0) {
			throw new ArithmeticException("Not an integer");
		}
		if (t.lowestPower > 0) {
			return reduced.zero();
		}
		return reduced.getElement(t.value);
	}

	public ModularIntegerRingElement reduce(PAdicNumber n, int numDigits) {
		PAdicNumber t = (PAdicNumber) n;
		if (t.lowestPower < 0) {
			throw new ArithmeticException("Not an integer");
		}
		ModularIntegerRing ring = reduction(numDigits);
		return ring.getElement(t.value.multiply(prime.pow(t.lowestPower)));
	}
	
	@Override
	public PAdicNumber upToUniformizerPower(PAdicNumber t) {
		return getInteger(t.value);
	}

	@Override
	public PAdicNumber lift(PFE t) {
		return new PAdicNumber(0, t.getValue());
	}

	public IntE roundToInteger(PAdicNumber t, int accuracy) {
		if (t.lowestPower < 0) {
			throw new ArithmeticException("Not an integer");
		}
		t = round(t, accuracy);
		BigInteger result = t.value.multiply(prime.pow(t.lowestPower));
		BigInteger biggestPower = prime.pow(accuracy);
		BigInteger alternateResult = result.subtract(biggestPower);
		return Integers.z().getInteger(result.abs().compareTo(alternateResult.abs()) < 0 ? result : alternateResult);
	}

	@Override
	public PAdicNumber round(PAdicNumber t, int accuracy) {
		if (t.lowestPower > accuracy) {
			return zero();
		}
		return new PAdicNumber(t.lowestPower, t.value.mod(prime.pow(accuracy - t.lowestPower)));
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<PAdicNumber> t) {
		return ringOfIntegers().isIrreducible(t);
	}

	@Override
	public FactorizationResult<Polynomial<PAdicNumber>> factorization(UnivariatePolynomial<PAdicNumber> t) {
		return ringOfIntegers().factorization(t, true, getAccuracy());
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public PAdicNumber add(PAdicNumber t1, PAdicNumber t2) {
		int lowestPower = Math.min(t1.lowestPower, t2.lowestPower);
		BigInteger result = t1.value.multiply(prime.pow(t1.lowestPower - lowestPower));
		result = result.add(t2.value.multiply(prime.pow(t2.lowestPower - lowestPower)));
		return new PAdicNumber(lowestPower, result);
	}

	@Override
	public PAdicNumber negative(PAdicNumber t) {
		return new PAdicNumber(t.lowestPower, prime.pow(accuracy - Math.min(t.lowestPower, 0)).subtract(t.value));
	}

	@Override
	public PAdicNumber negative(PAdicNumber t, int accuracy) {
		if (accuracy < this.accuracy) {
			return round(negative(t), accuracy);
		} else if (accuracy == this.accuracy) {
			return negative(t);
		}
		return withAccuracy(accuracy).negative(t);
	}

	@Override
	public PAdicNumber multiply(PAdicNumber t1, PAdicNumber t2) {
		return new PAdicNumber(t1.lowestPower + t2.lowestPower, t1.value.multiply(t2.value));
	}

	@Override
	public PAdicNumber inverse(PAdicNumber t) {
		return new PAdicNumber(-t.lowestPower, t.value.modInverse(prime.pow(accuracy + Math.abs(t.lowestPower))));
	}

	@Override
	public PAdicNumber inverse(PAdicNumber t, int accuracy) {
		if (accuracy < this.accuracy) {
			return round(inverse(t), accuracy);
		} else if (accuracy == this.accuracy) {
			return inverse(t);
		}
		return withAccuracy(accuracy).inverse(t);
	}

	@Override
	public PAdicNumber divide(PAdicNumber dividend, PAdicNumber divisor) {
		return new PAdicNumber(dividend.lowestPower - divisor.lowestPower, dividend.value.multiply(divisor.value
				.modInverse(prime.pow(accuracy + Math.abs(dividend.lowestPower) + Math.abs(divisor.lowestPower)))));
	}

	@Override
	public PivotStrategy<PAdicNumber> preferredPivotStrategy() {
		return new ValuationPivotStrategy<>(this.valuation());
	}

	public PAdicNumber rightShift(PAdicNumber t, int shift) {
		return new PAdicNumber(t.lowestPower - shift, t.value);
	}

	public PAdicNumber leftShift(PAdicNumber t, int shift) {
		return new PAdicNumber(t.lowestPower + shift, t.value);
	}

	@Override
	public PAdicNumber getInteger(BigInteger t) {
		boolean negative = t.compareTo(BigInteger.ZERO) < 0;
		if (negative) {
			t = t.negate();
		}
		PAdicNumber number = new PAdicNumber(0, t);
		return negative ? negative(number) : number;
	}

	@Override
	public String toString() {
		return "Q_" + prime;
	}

	@Override
	public PAdicNumber uniformizer() {
		return new PAdicNumber(1, BigInteger.ONE);
	}

	@Override
	public LocalRing<PAdicNumber, PFE> ringOfIntegers() {
		return localRing;
	}

	@Override
	public boolean isInteger(PAdicNumber t) {
		return valuation(t).compareTo(new Value(0)) >= 0;
	}
}
