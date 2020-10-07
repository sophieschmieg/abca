package fields.padics;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.finitefields.ModularIntegerRing;
import fields.finitefields.ModularIntegerRing.ModularIntegerRingElement;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PrimeFieldElement;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Field;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.ValueField;
import fields.padics.PAdicExtensionField.PAdicNumberExt;
import fields.padics.PAdicField.PAdicNumber;
import fields.polynomials.AbstractPolynomialRing;
import util.MiscAlgorithms;

public class PAdicField extends AbstractField<PAdicNumber> implements Field<PAdicNumber>, ValueField<PAdicNumber> {
	private BigInteger prime;
	private PrimeField reduced;
	private int maxAccuracy;
	final public static int INFINITE_ACCURACY = -1;

	public class PAdicNumber extends AbstractElement<PAdicNumber> {
		private List<BigInteger> digits;
		private List<BigInteger> periodDigits;
		private int lowestPower;
		private int accuracy;

		private PAdicNumber(int lowestPower, List<BigInteger> digits, List<BigInteger> periodDigits, int accuracy) {
			periodDigits = minimalPeriod(periodDigits);
			while (digits.size() > 0
					&& digits.get(digits.size() - 1).equals(periodDigits.get(periodDigits.size() - 1))) {
				try {
					digits.remove(digits.size() - 1);
				} catch (UnsupportedOperationException e) {
					List<BigInteger> newDigits = new ArrayList<>();
					newDigits.addAll(digits.subList(0, digits.size() - 1));
					digits = newDigits;
				}
				rotatePeriod(periodDigits);
			}
			while (digits.size() > 0 && digits.get(0).equals(BigInteger.ZERO)) {
				lowestPower++;
				try {
					digits.remove(0);
				} catch (UnsupportedOperationException e) {
					List<BigInteger> newDigits = new ArrayList<>();
					newDigits.addAll(digits.subList(1, digits.size()));
					digits = newDigits;
				}
				if (accuracy != INFINITE_ACCURACY) {
					accuracy--;
				}
			}
			while (digits.size() == 0 && periodDigits.size() > 1 && periodDigits.get(0).equals(BigInteger.ZERO)) {
				lowestPower++;
				rotatePeriod(periodDigits);
				if (accuracy != INFINITE_ACCURACY) {
					accuracy--;
				}
			}
			if (digits.size() == 0 && periodDigits.size() == 1 && periodDigits.get(0).equals(BigInteger.ZERO)) {
				lowestPower = 0;
			}
			this.lowestPower = lowestPower;
			this.periodDigits = periodDigits;
			this.digits = digits;
			this.accuracy = accuracy;
		}

		public String toString() {
			StringBuilder build = new StringBuilder().append("...[");
			for (int i = periodDigits.size() - 1; i >= 0; i--) {
				build.append(periodDigits.get(i).toString());
				if (prime.compareTo(BigInteger.TEN) > 0 && i != 0) {
					build.append(" ");
				}
			}
			build.append("]");
			for (int i = digits.size() - 1; i >= 0; i--) {
				build.append(digits.get(i).toString());
				if (prime.compareTo(BigInteger.TEN) > 0 && i != 0) {
					build.append(" ");
				}
			}
			build.append("*" + prime.toString() + "^" + lowestPower);
			return build.toString();
		}

		private BigInteger getDigit(int pos) {
			if (pos < lowestPower) {
				return BigInteger.ZERO;
			}
			if (pos < lowestPower + digits.size()) {
				return digits.get(pos - lowestPower);
			}
			return periodDigits.get((pos - lowestPower - digits.size()) % periodDigits.size());
		}

		@Override
		public int compareTo(PAdicNumber o) {
			boolean thisZero = false;
			if (this.digits.size() == 0 && this.periodDigits.size() == 1
					&& this.periodDigits.get(0).equals(BigInteger.ZERO)) {
				thisZero = true;
			}
			boolean otherZero = false;
			if (o.digits.size() == 0 && o.periodDigits.size() == 1 && o.periodDigits.get(0).equals(BigInteger.ZERO)) {
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
			if (this.lowestPower != o.lowestPower) {
				return this.lowestPower - o.lowestPower;
			}
			for (int i = lowestPower; i < lowestPower + Math.max(this.digits.size() + this.periodDigits.size(),
					o.digits.size() + o.periodDigits.size()); i++) {
				int cmp = this.getDigit(i).compareTo(o.getDigit(i));
				if (cmp != 0) {
					return cmp;
				}
			}
			return this.periodDigits.size() - o.periodDigits.size();
		}
	}

	private static List<BigInteger> minimalPeriod(List<BigInteger> period) {
		if (period.size() == 0) {
			ArrayList<BigInteger> result = new ArrayList<>();
			result.add(BigInteger.ZERO);
			return result;
		}
		periodLoop: for (int i = 1; i < period.size(); i++) {
			if (period.size() % i == 0) {
				for (int j = i; j < period.size(); j++) {
					if (!period.get(j).equals(period.get(j % i))) {
						continue periodLoop;
					}
				}
				return period.subList(0, i);
			}
		}
		return period;
	}

	private static void rotatePeriod(List<BigInteger> period) {
		BigInteger last = period.get(period.size() - 1);
		period.remove(period.size() - 1);
		period.add(0, last);
	}

	private static int minAccuracy(int accuracy1, int accuracy2) {
		if (accuracy1 == INFINITE_ACCURACY) {
			return accuracy2;
		}
		if (accuracy2 == INFINITE_ACCURACY) {
			return accuracy1;
		}
		return Math.min(accuracy1, accuracy2);
	}

	public PAdicField(BigInteger prime, int maxAccuracy) {
		if (!prime.isProbablePrime(100)) {
			throw new ArithmeticException("not a prime!");
		}
		this.maxAccuracy = maxAccuracy;
		this.prime = prime;
		this.reduced = new PrimeField(prime);
	}

	public int getMaxAccuracy() {
		return maxAccuracy;
	}

	public PrimeField reduction() {
		return reduced;
	}

	@Override
	public PAdicNumber getRandomElement() {
		Random r = new Random();
		List<BigInteger> digits = new ArrayList<>();
		List<BigInteger> period = new ArrayList<>();
		int lowestPower = r.nextInt(10) - 5;
		int digitsSize = r.nextInt(10);
		int periodSize = r.nextInt(5);
		for (int i = 0; i < digitsSize; i++) {
			digits.add(reduced.getRandomElement().getValue());
		}
		for (int i = 0; i < periodSize; i++) {
			period.add(reduced.getRandomElement().getValue());
		}
		return new PAdicNumber(lowestPower, digits, period, maxAccuracy);
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
		return new PAdicNumber(0, Collections.emptyList(), Collections.emptyList(), INFINITE_ACCURACY);
	}

	@Override
	public PAdicNumber one() {
		return new PAdicNumber(0, Collections.singletonList(BigInteger.ONE), Collections.emptyList(),
				INFINITE_ACCURACY);
	}

	@Override
	public double value(PAdicNumber t) {
		if (t.equals(zero())) {
			return 0.0;
		}
		double value = prime.pow(Math.abs(t.lowestPower)).doubleValue();
		if (t.lowestPower > 0) {
			return 1.0 / value;
		}
		return value;
	}

	public Fraction toRational(PAdicNumber t) {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		Fraction p = q.getEmbedding(z.getInteger(prime));
		Fraction nonPeriod = q.zero();
		int size = t.lowestPower + t.digits.size();
		for (int i = t.lowestPower; i < size; i++) {
			nonPeriod = q.add(nonPeriod, q.multiply(t.getDigit(i), q.power(p, i)));
		}
		Fraction numerator = q.zero();
		for (int i = size; i < size + t.periodDigits.size(); i++) {
			numerator = q.add(numerator, q.multiply(t.getDigit(i), q.power(p, i)));
		}
		Fraction denominator = q.subtract(q.one(), q.power(p, t.periodDigits.size()));
		Fraction period = q.divide(numerator, denominator);
		return q.add(nonPeriod, period);
	}

	public Polynomial<Fraction> toRational(Polynomial<PAdicNumber> t) {
		PolynomialRing<PAdicNumber> ring = t.getPolynomialRing();
		PolynomialRing<Fraction> rationalRing = AbstractPolynomialRing.getPolynomialRing(Rationals.q(),
				ring.numberOfVariables(), ring.getComparator());
		return MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public Fraction evaluate(PAdicNumber n) {
				return toRational(n);
			}
		}, rationalRing);
	}

	public PAdicNumber fromRational(Fraction t) {
		BigInteger numerator = t.getNumerator().getValue();
		BigInteger denominator = t.getDenominator().getValue();
		return divide(getInteger(numerator), getInteger(denominator));
	}

	public Polynomial<PAdicNumber> fromRational(Polynomial<Fraction> t) {
		PolynomialRing<Fraction> ring = t.getPolynomialRing();
		PolynomialRing<PAdicNumber> padicRing = AbstractPolynomialRing.getPolynomialRing(this, ring.numberOfVariables(),
				ring.getComparator());
		return MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public PAdicNumber evaluate(Fraction n) {
				return fromRational(n);
			}
		}, padicRing);
	}

	public PrimeFieldElement reduce(PAdicNumber t) {
		if (t.lowestPower < 0) {
			throw new ArithmeticException("Not an integer");
		}
		return reduced.getElement(t.getDigit(0));
	}

	public Polynomial<PrimeFieldElement> reduce(Polynomial<PAdicNumber> t) {
		PolynomialRing<PAdicNumber> ring = t.getPolynomialRing();
		PolynomialRing<PrimeFieldElement> reducedRing = AbstractPolynomialRing.getPolynomialRing(reduced,
				ring.numberOfVariables(), ring.getComparator());
		return MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public PrimeFieldElement evaluate(PAdicNumber n) {
				return reduce(n);
			}
		}, reducedRing);
	}

	public ModularIntegerRingElement reduce(PAdicNumber t, int numDigits) {
		if (t.lowestPower < 0) {
			throw new ArithmeticException("Not an integer");
		}
		ModularIntegerRing ring = new ModularIntegerRing(prime.pow(numDigits));
		ModularIntegerRingElement result = ring.zero();
		ModularIntegerRingElement p = ring.getElement(prime);
		for (int i = 0; i < numDigits; i++) {
			result = ring.add(result, ring.multiply(t.getDigit(i), ring.power(p, i)));
		}
		return result;
	}

	public PAdicNumber lift(PrimeFieldElement t) {
		return new PAdicNumber(0, Collections.singletonList(t.getValue()), Collections.emptyList(), INFINITE_ACCURACY);
	}

	public Polynomial<PAdicNumber> lift(Polynomial<PrimeFieldElement> t) {
		PolynomialRing<PrimeFieldElement> ring = t.getPolynomialRing();
		PolynomialRing<PAdicNumber> padicRing = AbstractPolynomialRing.getPolynomialRing(this, ring.numberOfVariables(),
				ring.getComparator());
		return MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public PAdicNumber evaluate(PrimeFieldElement n) {
				return lift(n);
			}
		}, padicRing);
	}

	public IntE roundToInteger(PAdicNumber t, int accuracy) {
		if (t.lowestPower < 0) {
			throw new ArithmeticException("Not an integer");
		}
		BigInteger result = BigInteger.ZERO;
		for (int i = 0; i < accuracy; i++) {
			result = result.add(t.getDigit(i).multiply(prime.pow(i)));
		}
		BigInteger biggestPower = prime.pow(accuracy);
		BigInteger alternateResult = result.subtract(biggestPower);
		return Integers.z().getInteger(result.abs().compareTo(alternateResult.abs()) < 0 ? result : alternateResult);
	}

	public PAdicNumber round(PAdicNumber t, int accuracy) {
		if (accuracy == INFINITE_ACCURACY) {
			return t;
		}
		List<BigInteger> digits = new ArrayList<>();
		for (int i = t.lowestPower; i < accuracy; i++) {
			digits.add(t.getDigit(i));
		}
		return new PAdicNumber(t.lowestPower, digits, Collections.emptyList(), INFINITE_ACCURACY);
	}

	public Polynomial<PAdicNumber> round(Polynomial<PAdicNumber> t, int accuracy) {
		return MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public PAdicNumber evaluate(PAdicNumber n) {
				return round(n, accuracy);
			}
		}, t.getPolynomialRing());
	}

	@Override
	public List<Polynomial<PAdicNumber>> factorization(Polynomial<PAdicNumber> t) {
		PAdicExtensionField base = new PAdicExtensionField(this);
		List<Polynomial<PAdicNumberExt>> factors = base
				.factorization(MiscAlgorithms.mapPolynomial(t, new MathMap<PAdicNumber, PAdicNumberExt>() {
					@Override
					public PAdicNumberExt evaluate(PAdicNumber t) {
						return base.getEmbedding(t);
					}
				}, base.getUnivariatePolynomialRing()));
		List<Polynomial<PAdicNumber>> result = new ArrayList<>();
		for (Polynomial<PAdicNumberExt> factor : factors) {
			result.add(MiscAlgorithms.mapPolynomial(factor, new MathMap<PAdicNumberExt, PAdicNumber>() {
				@Override
				public PAdicNumber evaluate(PAdicNumberExt t) {
					return t.asExtensionFieldElement().getElement().get(0);
				}
			}, this.getUnivariatePolynomialRing()));
		}
		return result;
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public PAdicNumber add(PAdicNumber t1, PAdicNumber t2) {
		int lowestPower = Math.min(t1.lowestPower, t2.lowestPower);
		int size = Math.max(t1.lowestPower + t1.digits.size(), t2.lowestPower + t2.digits.size())
				+ Math.max(2 * t1.periodDigits.size(), 2 * t2.periodDigits.size());
		int accuracy = minAccuracy(t1.accuracy, t2.accuracy);
		List<BigInteger> digits = new ArrayList<>();
		List<BigInteger> period = new ArrayList<>();
		BigInteger carry = BigInteger.ZERO;
		for (int i = lowestPower; i < size; i++) {
			BigInteger result = t1.getDigit(i).add(t2.getDigit(i)).add(carry);
			digits.add(result.mod(prime));
			carry = result.subtract(result.mod(prime)).divide(prime);
		}
		BigInteger initialCarry = carry;
		int periodSize = t1.periodDigits.size() * t2.periodDigits.size();
		for (int i = size; i < size + periodSize; i++) {
			BigInteger result = t1.getDigit(i).add(t2.getDigit(i)).add(carry);
			period.add(result.mod(prime));
			carry = result.subtract(result.mod(prime)).divide(prime);
		}
		if (!initialCarry.equals(carry)) {
			throw new ArithmeticException("Algorithm wrong");
		}
		PAdicNumber result = new PAdicNumber(lowestPower, digits, period, accuracy);
		if (this.maxAccuracy != INFINITE_ACCURACY
				&& result.digits.size() + result.periodDigits.size() > this.maxAccuracy) {
			result = this.round(result, this.maxAccuracy);
		}
		return result;
	}

	@Override
	public PAdicNumber negative(PAdicNumber t) {
		BigInteger one = BigInteger.ONE;
		BigInteger pmone = prime.subtract(one);
		List<BigInteger> digits = new ArrayList<>();
		List<BigInteger> period = new ArrayList<>();
		for (int i = 0; i < t.digits.size(); i++) {
			digits.add(pmone.subtract(t.digits.get(i)));
		}
		for (int i = 0; i < t.periodDigits.size(); i++) {
			period.add(pmone.subtract(t.periodDigits.get(i)));
		}
		return add(new PAdicNumber(t.lowestPower, Collections.singletonList(one), Collections.emptyList(), t.accuracy),
				new PAdicNumber(t.lowestPower, digits, period, t.accuracy));
	}

	@Override
	public PAdicNumber multiply(PAdicNumber t1, PAdicNumber t2) {
		if (t2.equals(zero())) {
			return zero();
		}
		if (maxAccuracy == INFINITE_ACCURACY) {
			return divide(t1, inverse(t2));
		}
		PAdicNumber result = zero();
		for (int i = t1.lowestPower; i <= maxAccuracy; i++) {
			result = add(result, leftShift(multiplyDigit(t1.getDigit(i), t2), i));
		}
		return round(result, maxAccuracy);
	}

	@Override
	public PAdicNumber inverse(PAdicNumber t) {
		return divide(one(), t);
	}

	@Override
	public PAdicNumber divide(PAdicNumber dividend, PAdicNumber divisor) {
		List<BigInteger> digits = new ArrayList<>();
		Map<PAdicNumber, Integer> residue = new TreeMap<>();
		int dividendPower = dividend.lowestPower;
		PAdicNumber alignedDividend = rightShift(dividend, dividendPower);
		int divisorPower = divisor.lowestPower;
		PAdicNumber alignedDivisor = rightShift(divisor, divisorPower);
		int accuracy = minAccuracy(dividend.accuracy, divisor.accuracy);
		BigInteger divisorMultiplier = reduced.inverse(reduce(alignedDivisor)).getValue();
		for (int i = 0; !residue.containsKey(alignedDividend); i++) {
			residue.put(alignedDividend, i);
			if ((maxAccuracy != INFINITE_ACCURACY && i > maxAccuracy)
					|| (accuracy != INFINITE_ACCURACY && i > accuracy)) {
				accuracy = minAccuracy(accuracy, maxAccuracy);
				break;
			}
			PAdicNumber result = multiplyDigit(divisorMultiplier, alignedDividend);
			BigInteger digit = result.getDigit(0);
			digits.add(digit);
			alignedDividend = rightShift(subtract(alignedDividend, multiplyDigit(digit, alignedDivisor)), 1);
		}
		int periodIndex = residue.get(alignedDividend);
		List<BigInteger> nonPeriod = new ArrayList<>(digits.subList(0, periodIndex));
		List<BigInteger> period = new ArrayList<>(digits.subList(periodIndex, digits.size()));
		PAdicNumber result = new PAdicNumber(dividendPower - divisorPower, nonPeriod, period, accuracy);
		if (this.maxAccuracy != INFINITE_ACCURACY
				&& result.digits.size() + result.periodDigits.size() > this.maxAccuracy) {
			result = this.round(result, this.maxAccuracy);
		}
		return result;
	}

	public PAdicNumber rightShift(PAdicNumber t, int shift) {
		return new PAdicNumber(t.lowestPower - shift, t.digits, t.periodDigits, t.accuracy);
	}

	public PAdicNumber leftShift(PAdicNumber t, int shift) {
		return new PAdicNumber(t.lowestPower + shift, t.digits, t.periodDigits, t.accuracy);
	}

	private PAdicNumber multiplyDigit(BigInteger digit, PAdicNumber t) {
		if (digit.equals(BigInteger.ZERO)) {
			return zero();
		}
		List<BigInteger> digits = new ArrayList<>();
		List<BigInteger> period = new ArrayList<>();
		BigInteger carry = BigInteger.ZERO;
		int size = t.lowestPower + t.digits.size() + 2 * t.periodDigits.size();
		for (int i = t.lowestPower; i < size; i++) {
			BigInteger result = t.getDigit(i).multiply(digit).add(carry);
			digits.add(result.mod(prime));
			carry = result.subtract(result.mod(prime)).divide(prime);
		}
		BigInteger initialCarry = carry;
		for (int i = size; i < size + t.periodDigits.size(); i++) {
			BigInteger result = t.getDigit(i).multiply(digit).add(carry);
			period.add(result.mod(prime));
			carry = result.subtract(result.mod(prime)).divide(prime);
		}
		if (!carry.equals(initialCarry)) {
			throw new ArithmeticException("algorithm wrong");
		}
		return new PAdicNumber(t.lowestPower, digits, period, t.accuracy);
	}

	@Override
	public PAdicNumber getInteger(BigInteger t) {
		boolean negative = t.compareTo(BigInteger.ZERO) < 0;
		if (negative) {
			t = t.negate();
		}
		List<BigInteger> digits = new ArrayList<>();
		while (!t.equals(BigInteger.ZERO)) {
			digits.add(t.mod(prime));
			t = t.subtract(t.mod(prime)).divide(prime);
		}
		PAdicNumber number = new PAdicNumber(0, digits, Collections.emptyList(), INFINITE_ACCURACY);
		return negative ? negative(number) : number;
	}

	@Override
	public String toString() {
		return "Z_" + prime;
	}
}
