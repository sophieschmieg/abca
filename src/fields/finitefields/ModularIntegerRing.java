package fields.finitefields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

import fields.exceptions.InfinityException;
import fields.finitefields.ModularIntegerRing.ModularIntegerRingElement;
import fields.helper.AbstractIdeal;
import fields.helper.AbstractRing;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
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
		return this.getElement(t1.value.add(t2.value));
	}

	@Override
	public ModularIntegerRingElement negative(ModularIntegerRingElement t) {
		return this.getElement(n.subtract(t.value));
	}

	@Override
	public ModularIntegerRingElement multiply(ModularIntegerRingElement t1, ModularIntegerRingElement t2) {
		return this.getElement(t1.value.multiply(t2.value));
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
		Map<BigInteger, Integer> primes = MiscAlgorithms.primeDecomposition(n);
		BigInteger result = BigInteger.ONE;
		for (BigInteger p : primes.keySet()) {
			result = result.multiply(p.pow(primes.get(p) - 1));
			result = result.multiply(p.subtract(BigInteger.ONE));
		}
		return result;
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

	@Override
	public boolean isDivisible(ModularIntegerRingElement dividend, ModularIntegerRingElement divisor) {
		if (isUnit(divisor)) {
			return true;
		}
		return isZeroDivisor(dividend);
	}

	@Override
	public List<ModularIntegerRingElement> quotientAndRemainder(ModularIntegerRingElement dividend,
			ModularIntegerRingElement divisor) {
		if (!isUnit(divisor)) {
			throw new ArithmeticException("not divisble");
		}
		List<ModularIntegerRingElement> result = new ArrayList<>();
		result.add(multiply(dividend, inverse(divisor)));
		result.add(zero());
		return result;
	}

	@Override
	public BigInteger euclidMeasure(ModularIntegerRingElement t) {
		return null;
	}
	
	@Override
	public ModularIntegerRingElement projectToUnit(ModularIntegerRingElement t) {
		if (t.equals(zero())) {
			return one();
		}
		return getElement(t.value.divide(t.value.gcd(n)));
	}
	
	public Ideal<ModularIntegerRingElement> getIdeal(List<ModularIntegerRingElement> generators) {
		if (generators.size() == 0) {
			return new ModularIntegerIdeal(0);
		}
		Integers z = Integers.z();
		IntE generator = z.lift(generators.get(0));
		for (int i = 1; i < generators.size(); i++) {
			generator = z.gcd(generator, z.lift(generators.get(i)));
		}
		generator = z.gcd(generator, z.getInteger(n));
		return new ModularIntegerIdeal(generator.getValue());

	}
	
	public Ideal<ModularIntegerRingElement> intersect(Ideal<ModularIntegerRingElement> t1, Ideal<ModularIntegerRingElement> t2) {
		return multiply(t1, t2);
	}
	
	public int krullDimension() {
		return n.equals(BigInteger.ONE) ? -1 : 0;
	}

	public class ModularIntegerRingElement implements Element<ModularIntegerRingElement> {
		private BigInteger value;

		private ModularIntegerRingElement(BigInteger value) {
			this.value = value.mod(n);
			while (this.value.compareTo(BigInteger.ZERO) < 0)
				this.value = this.value.add(n);
		}

		public BigInteger getValue() {
			return value;
		}

		public boolean equals(Object o) {
			if (!(o instanceof ModularIntegerRingElement))
				return false;
			ModularIntegerRingElement e = (ModularIntegerRingElement) o;
			return this.value.equals(e.value);
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
			return getElement(t.value.mod(m.value));
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
	}
}
