package fields.finitefields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.finitefields.ModularIntegerRing.ModularIntegerRingElement;
import fields.helper.AbstractElement;
import fields.helper.AbstractIdeal;
import fields.helper.AbstractRing;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Ideal;
import fields.local.Value;
import util.MiscAlgorithms;

public class ModularIntegerRing extends AbstractRing<ModularIntegerRingElement> {
	private BigInteger n;

	public ModularIntegerRing(int n) {
		this.n = BigInteger.valueOf(n);
	}

	public ModularIntegerRing(BigInteger n) {
		this.n = n;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public ModularIntegerRingElement zero() {
		return this.getElement(0);
	}

	@Override
	public ModularIntegerRingElement one() {
		return this.getElement(1);
	}

	@Override
	public BigInteger characteristic() {
		return n;
	}

	@Override
	public ModularIntegerRingElement add(ModularIntegerRingElement t1, ModularIntegerRingElement t2) {
		return getElement(t1.value.add(t2.value));
	}

	@Override
	public ModularIntegerRingElement negative(ModularIntegerRingElement t) {
		return this.getElement(n.subtract(t.value));
	}

	@Override
	public ModularIntegerRingElement multiply(ModularIntegerRingElement t1, ModularIntegerRingElement t2) {
		return getElement(t1.value.multiply(t2.value));
	}

	@Override
	public ModularIntegerRingElement inverse(ModularIntegerRingElement t) {
		return getElement(t.value.modInverse(n));
	}

	@Override
	public ModularIntegerRingElement getRandomElement() {
		return getElement(MiscAlgorithms.randomBigInteger(new Random(), n));
	}

	@Override
	public BigInteger getNumberOfElements() {
		return this.n;
	}

	public ModularIntegerRingElement reduce(IntE t) {
		return new ModularIntegerRingElement(t.getValue());
	}

	@Override
	public Iterator<ModularIntegerRingElement> iterator() {
		return new Iterator<ModularIntegerRingElement>() {
			private BigInteger i = BigInteger.ZERO;

			@Override
			public boolean hasNext() {
				return this.i.compareTo(n) < 0;
			}

			@Override
			public ModularIntegerRingElement next() {
				this.i = this.i.add(BigInteger.ONE);
				return ModularIntegerRing.this.getElement(i.subtract(BigInteger.ONE));
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	@Override
	public boolean isCommutative() {
		return true;
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public String toString() {
		return "Z/" + n + "Z";
	}

	public ModularIntegerRingElement getElement(int value) {
		return new ModularIntegerRingElement(BigInteger.valueOf(value));
	}

	public ModularIntegerRingElement getElement(BigInteger value) {
		return new ModularIntegerRingElement(value);
	}

	@Override
	public boolean isUnit(ModularIntegerRingElement t) {
		return n.gcd(t.value).equals(BigInteger.ONE);
	}

	@Override
	public BigInteger getNumberOfUnits() {
		Integers z = Integers.z();
		return z.eulerToitent(new IntE(n)).getValue();
	}

	@Override
	public Iterable<ModularIntegerRingElement> getUnits() {
		return new Iterable<ModularIntegerRingElement>() {

			@Override
			public Iterator<ModularIntegerRingElement> iterator() {
				return new Iterator<ModularIntegerRingElement>() {
					private BigInteger i = BigInteger.ZERO;

					@Override
					public boolean hasNext() {
						return this.i.compareTo(n) < 0;
					}

					@Override
					public ModularIntegerRingElement next() {
						do {
							this.i = this.i.add(BigInteger.ONE);
						} while (!ModularIntegerRing.this.isUnit(ModularIntegerRing.this.getElement(i)));
						return ModularIntegerRing.this.getElement(i);
					}
				};
			}
		};
	}

	@Override
	public boolean isIntegral() {
		return n.isProbablePrime(10);
	}

	@Override
	public boolean isZeroDivisor(ModularIntegerRingElement t) {
		return !isUnit(t);
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	public boolean isDedekindDomain() {
		return false;
	}

	@Override
	public boolean isDivisible(ModularIntegerRingElement dividend, ModularIntegerRingElement divisor) {
		BigInteger divisorGcd = divisor.value.gcd(n);
		BigInteger dividendGcd = dividend.value.gcd(n);
		return dividendGcd.mod(divisorGcd).equals(BigInteger.ZERO);
	}

	@Override
	public QuotientAndRemainderResult<ModularIntegerRingElement> quotientAndRemainder(
			ModularIntegerRingElement dividend, ModularIntegerRingElement divisor) {
		BigInteger divisorGcd = divisor.value.gcd(n);
		BigInteger dividendGcd = dividend.value.gcd(n);
		BigInteger[] qr = dividendGcd.divideAndRemainder(divisorGcd);
		if (!qr[1].equals(BigInteger.ZERO)) {
			throw new ArithmeticException("not divisble");
		}
		BigInteger reducedDividend = dividend.value.divide(divisorGcd);
		BigInteger reducedDivisor = divisor.value.divide(divisorGcd);
		BigInteger reducedN = n.divide(divisorGcd);
		ModularIntegerRingElement result = new ModularIntegerRingElement(reducedDividend.multiply(reducedDivisor.modInverse(reducedN)));
		return new QuotientAndRemainderResult<>(result, zero());
	}

	@Override
	public BigInteger euclidMeasure(ModularIntegerRingElement t) {
		return null;
	}

	@Override
	public FactorizationResult<ModularIntegerRingElement,ModularIntegerRingElement> uniqueFactorization(ModularIntegerRingElement t) {
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(z.getInteger(n));
		IntE e = z.getInteger(t.value);
		SortedMap<ModularIntegerRingElement, Integer> result = new TreeMap<>();
		for (IntE prime : factors.primeFactors()) {
			int power = factors.multiplicity(prime);
			IntE factor = z.power(prime, power);
			IntE gcd = z.gcd(factor, e);
			if (!z.isUnit(gcd)) {
				e = z.divideChecked(e, gcd);
				int factorPower = 0;
				while (!z.isUnit(gcd)) {
					factorPower++;
					gcd = z.divideChecked(gcd, prime);
				}
				result.put(getInteger(prime), factorPower);
			}
		}
		if (!isUnit(getInteger(e))) {
			throw new ArithmeticException("Algorithm wrong");
		}
		return new FactorizationResult<>(getInteger(e), result);
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return true;
	}

	@Override
	public ModularIntegerRingElement projectToUnit(ModularIntegerRingElement t) {
		if (t.equals(zero())) {
			return one();
		}
		return getElement(t.value.divide(t.value.gcd(n)));
	}

	public IdealResult<ModularIntegerRingElement, ModularIntegerIdeal> getIdealWithTransforms(List<ModularIntegerRingElement> generators) {
		if (generators.size() == 0) {
			return new IdealResult<>(Collections.singletonList(Collections.emptyList()), generators,
					new ModularIntegerIdeal(0));
		}
		Integers z = Integers.z();
		List<IntE> lifted = new ArrayList<>();
		for (ModularIntegerRingElement g : generators) {
			lifted.add(z.lift(g));
		}
		ExtendedEuclideanListResult<IntE> extendedEuclidean = z.extendedEuclidean(lifted);
		List<ModularIntegerRingElement> reduced = new ArrayList<>();
		for (IntE expression : extendedEuclidean.getCoeffs()) {
			reduced.add(reduce(expression));
		}

		return new IdealResult<>(Collections.singletonList(reduced), generators,
				new ModularIntegerIdeal(z.gcd(extendedEuclidean.getGcd(), z.getInteger(n)).getValue()));
	}

	public Ideal<ModularIntegerRingElement> intersect(Ideal<ModularIntegerRingElement> t1,
			Ideal<ModularIntegerRingElement> t2) {
		return getIdeal(Collections.singletonList(lcm(t1.generators().get(0), t2.generators().get(0))));
	}

	@Override
	public Ideal<ModularIntegerRingElement> radical(Ideal<ModularIntegerRingElement> t) {
		Integers z = Integers.z();
		IntE m = z.lift(t.generators().get(0));
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(m);
		ModularIntegerRingElement radical = one();
		for (IntE prime : factors.primeFactors()) {
			radical = multiply(radical, reduce(prime));
		}
		return getIdeal(Collections.singletonList(radical));
	}

	public int krullDimension() {
		return n.equals(BigInteger.ONE) ? -1 : 0;
	}

	public class ModularIntegerRingElement extends AbstractElement<ModularIntegerRingElement> {
		private BigInteger value;

		private ModularIntegerRingElement(BigInteger value) {
			this.value = value.mod(n);
			while (this.value.compareTo(BigInteger.ZERO) < 0)
				this.value = this.value.add(n);
		}

		public BigInteger getValue() {
			return value;
		}

		public String toString() {
			return this.value.toString();
		}

		@Override
		public int compareTo(ModularIntegerRingElement o) {
			return this.value.compareTo(o.value);
		}

	}

	public class ModularIntegerIdeal extends AbstractIdeal<ModularIntegerRingElement> {
		private ModularIntegerRingElement m;

		private ModularIntegerIdeal(ModularIntegerRingElement m) {
			super(ModularIntegerRing.this);
			this.m = m;
		}

		private ModularIntegerIdeal(int m) {
			this(BigInteger.valueOf(m));
		}

		private ModularIntegerIdeal(BigInteger m) {
			this(new ModularIntegerRingElement(m));
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
		public List<ModularIntegerRingElement> generators() {
			return Collections.singletonList(m);
		}

		@Override
		public List<ModularIntegerRingElement> generate(ModularIntegerRingElement t) {
			return Collections.singletonList(getElement(t.value.divide(m.value)));
		}

		@Override
		public ModularIntegerRingElement residue(ModularIntegerRingElement t) {
			return new ModularIntegerRingElement(t.value.mod(m.value));
		}

		@Override
		public boolean isPrime() {
			return m.value.isProbablePrime(100);
		}

		@Override
		public boolean isMaximal() {
			return isPrime() && !m.equals(ModularIntegerRing.this.zero());
		}

		@Override
		public boolean contains(ModularIntegerRingElement t) {
			return t.value.mod(m.value).equals(BigInteger.ZERO);
		}

		@Override
		public Value maximumPowerContains(ModularIntegerRingElement t) {
			if (t.equals(zero()) || m.equals(one())) {
				return Value.INFINITY;
			}
			if (m.equals(zero())) {
				return Value.ZERO;
			}
			FactorizationResult<Ideal<ModularIntegerRingElement>, Ideal<ModularIntegerRingElement>> factors = idealFactorization(this);
			Value value = Value.INFINITY;
			Integers z = Integers.z();
			for (Ideal<ModularIntegerRingElement> factor : factors.primeFactors()) {
				Ideal<IntE> liftedFactor = z.getIdeal(Collections.singletonList(z.lift(factor.generators().get(0))));
				value = value.min(new Value(
						z.valuation(z.lift(t), liftedFactor).value() / z.valuation(z.lift(m), liftedFactor).value()));
			}
			return value;
		}
	}
}