package fields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


public abstract class AbstractRing<T extends Element> implements Ring<T> {
	public Group<T> getAdditiveGroup() {
		return new Group<T>() {

			@Override
			public T getRandomElement() {
				return AbstractRing.this.getRandomElement();
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				return AbstractRing.this.getNumberOfElements();
			}

			@Override
			public Iterator<T> iterator() {
				return AbstractRing.this.iterator();
			}

			@Override
			public T neutral() {
				return AbstractRing.this.zero();
			}

			@Override
			public T inverse(T t) {
				return AbstractRing.this.negative(t);
			}

			@Override
			public T operate(T t1, T t2) {
				return AbstractRing.this.add(t1, t2);
			}

			@Override
			public boolean isFinite() {
				return AbstractRing.this.isFinite();
			}

		};
	}
	public Group<T> getMultiplicativeGroup() {
		return new Group<T>() {

			@Override
			public T getRandomElement() {
				T t;
				do {
					t = AbstractRing.this.getRandomElement();
				} while (!AbstractRing.this.isUnit(t));
				return t;
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				return AbstractRing.this.getNumberOfUnits();
			}

			@Override
			public Iterator<T> iterator() {
				return AbstractRing.this.getUnits().iterator();
			}

			@Override
			public T neutral() {
				return AbstractRing.this.one();
			}

			@Override
			public T inverse(T t) {
				return AbstractRing.this.inverse(t);
			}

			@Override
			public T operate(T t1, T t2) {
				return AbstractRing.this.multiply(t1, t2);
			}

			@Override
			public boolean isFinite() {
				return AbstractRing.this.isFinite();
			}
		};
	}

	public T add(T t1, T t2, T t3) {
		return this.add(this.add(t1, t2), t3);
	}

	public T subtract(T minuend, T subtrahend) {
		return this.add(minuend, this.negative(subtrahend));
	}
	public T getInteger(int n) {
		return multiply(n, this.one());
	}
	@Override
	public T getInteger(BigInteger n) {
		return multiply(n, this.one());
	}
	@Override
	public T multiply(int n, T t) {
		return this.multiply(BigInteger.valueOf(n), t);
	}
	@Override
	public T multiply(BigInteger n, T t) {
		T result = this.zero();
		if (n.signum() < 0) {
			n = n.negate();
			t = this.negative(t);
		}
		for (int i = 0; i <= n.bitLength(); i++) {
			if (n.testBit(i))
				result = this.add(result, t);
			t = this.add(t, t);
		}
		return result;
	}
	@Override
	public T multiply(T t1, T t2, T t3) {
		return this.multiply(t1, this.multiply(t2, t3));
	}
	@Override
	public T multiply(int n, T t1, T t2) {
		return this.multiply(n, this.multiply(t1, t2));
	}
	@Override
	public T multiply(BigInteger n, T t1, T t2) {
		return this.multiply(n, this.multiply(t1, t2));
	}
	@Override
	public T multiply(int n, T t1, T t2, T t3) {
		return this.multiply(n, this.multiply(t1, t2, t3));
	}
	@Override
	public T multiply(BigInteger n, T t1, T t2, T t3) {
		return this.multiply(n, this.multiply(t1, t2, t3));
	}

	@Override
	public T power(T t, int n) {
		return this.power(t, BigInteger.valueOf(n));
	}
	@Override
	public T power(T t, BigInteger n) {
		if (n.equals(BigInteger.ONE)) {
			return t;
		} else if (n.equals(BigInteger.TWO)) {
			return this.multiply(t, t);
		}
		T result;
		if (n.signum() < 0) {
			n = n.negate();
			t = this.inverse(t);
		}
		if (n.testBit(0)) {
			result = t;
		} else {
			result = this.one();
		}
		for (int i = 1; i <= n.bitLength(); i++) {
			t = this.multiply(t, t);
			if (n.testBit(i)) {
				result = this.multiply(result, t);
			}
		}
		return result;
	}

	public T gcd(T t1, T t2) {
		if (!this.isEuclidean())
			throw new RuntimeException("Not possible");
		return this.extendedEuclidean(t1, t2).get(0);
	}
	public List<T> extendedEuclidean(T t1, T t2) {
		if (!this.isEuclidean())
			throw new RuntimeException("Not possible");
		if (t1.equals(this.zero()) && t2.equals(this.zero())) {
			List<T> result = new ArrayList<T>();
			result.add(this.zero());
			result.add(this.zero());
			result.add(this.zero());
			return result;
		}
		if (t1.equals(this.zero())) {
			List<T> result = new ArrayList<T>();
			result.add(t2);
			result.add(this.zero());
			result.add(this.one());
			return result;
		}
		if (t2.equals(this.zero())) {
			List<T> result = new ArrayList<T>();
			result.add(t1);
			result.add(this.one());
			result.add(this.zero());
			return result;
		}
		T u, v, xp, x, yp, y, q, r;
		u = t1;
		v = t2;
		xp = this.one();
		yp = this.zero();
		x = this.zero();
		y = this.one();
		// x(0) * t1 + y(0) * t2 = u = r(0)
		// x(1) * t1 + y(1) * t2 = v = r(1)
		while (true) {
			List<T> qr = this.quotientAndRemainder(u, v);
			q = qr.get(0);
			r = qr.get(1);
			if (r.equals(this.zero()))
				break;
			// u = q * v + r
			// r(n+1) = r(n-1) - q * r(n)
			// r(n+1) = (x(n-1) - q * x(n)) * t1 + (y(n-1) - q * y(n)) * t2
			// x(n+1) = x(n-1) - q * x(n)
			// y(n+1) = y(n-1) - q * y(n)
			u = v;
			v = r;
			r = x;
			x = this.subtract(xp, this.multiply(q, x));
			xp = r;
			r = y;
			y = this.subtract(yp, this.multiply(q, y));
			yp = r;
		}
		List<T> result = new ArrayList<T>();
		result.add(v);
		result.add(x);
		result.add(y);
		return result;
	}

	@Override
	public BigInteger getNumberOfUnits() {
		return this.getNumberOfUnits();
	}
}
