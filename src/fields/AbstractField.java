package fields;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public abstract class AbstractField<T extends Element> implements Field<T>, Ring<T> {
	private T primitiveRoot = null;

	@Override
	public Group<T> getAdditiveGroup() {
		return new Group<T>() {

			@Override
			public T getRandomElement() {
				return this.getRandomElement();
			}

			@Override
			public int getNumberOfElements() throws InfinityException {
				return this.getNumberOfElements();
			}

			@Override
			public Iterable<T> getElements() throws InfinityException {
				return this.getElements();
			}

			@Override
			public T neutral() {
				return AbstractField.this.zero();
			}

			@Override
			public T inverse(T t) {
				return AbstractField.this.negative(t);
			}

			@Override
			public T operate(T t1, T t2) {
				return AbstractField.this.add(t1, t2);
			}

			@Override
			public boolean isFinite() {
				return this.isFinite();
			}
			
		};
	}

	@Override
	public Group<T> getMultiplicativeGroup() {
		return new Group<T>() {

			@Override
			public T getRandomElement() {
				T t;
				do {
					t = this.getRandomElement();
				} while (t.equals(AbstractField.this.zero()));
				return t;
			}

			@Override
			public int getNumberOfElements() throws InfinityException {
				return this.getNumberOfElements() - 1;
			}

			@Override
			public Iterable<T> getElements() throws InfinityException {
				return AbstractField.this.getNonZeroElements();
			}

			@Override
			public T neutral() {
				return AbstractField.this.one();
			}

			@Override
			public T inverse(T t) {
				return AbstractField.this.inverse(t);
			}

			@Override
			public T operate(T t1, T t2) {
				return AbstractField.this.multiply(t1, t2);
			}

			@Override
			public boolean isFinite() {
				return AbstractField.this.isFinite();
			}
		};
	}

	@Override
	public T add(T t1, T t2, T t3) {
		return this.add(this.add(t1, t2), t3);
	}

	@Override
	public T subtract(T minuend, T subtrahend) {
		return this.add(minuend, this.negative(subtrahend));
	}
	@Override
	public T getInteger(int n) {
		return multiply(n, this.one());
	}
	@Override
	public T multiply(int n, T t) {
		T result = this.zero();
		if (n < 0) {
			n = -n;
			t = this.negative(t);
		}
		while (n > 0) {
			if ((n & 1) == 1)
				result = this.add(result, t);
			n >>= 1;
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
	public T multiply(int n, T t1, T t2, T t3) {
		return this.multiply(n, this.multiply(t1, t2, t3));
	}
	
	@Override
	public T divide(T dividend, T divisor) {
		return this.multiply(dividend, this.inverse(divisor));
	}
	
	@Override
	public T power(T t, int n) {
		T result = this.one();
		if (n < 0) {
			n = -n;
			t = this.inverse(t);
		}
		while (n > 0) {
			if ((n & 1) == 1)
				result = this.multiply(result, t);
			n >>= 1;
			t = this.multiply(t, t);
		}
		return result;
	}

	public Iterable<T> getNonZeroElements() throws InfinityException {
		return new Iterable<T>() {
			private Iterable<T> elements = AbstractField.this.getElements();
			@Override
			public Iterator<T> iterator() {
				return new Iterator<T>() {
					private Iterator<T> it = elements.iterator();

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public T next() {
						T result = null;
						do {
							if (it.hasNext())
								result = it.next();
						} while (result.equals(AbstractField.this.zero()));
						return result;
					}

					@Override
					public void remove() {
						it.remove();
					}
				};
			}
			
		};
	}

	@Override
	public boolean isUnit(T t) {
		return !t.equals(this.zero());
	}

	@Override
	public boolean hasRoot(T t, int n) {
		if (t.equals(this.zero()))
			return true;
		if (!this.isFinite())
			throw new RuntimeException("Not implemented");
		int q = this.getNumberOfElements() - 1;
		int u = q;
		int v = n;
		while (u % v != 0) {
			int r = u % v;
			u = v;
			v = r;
		}
		T p = this.power(t, q / v);
		return p.equals(this.one());
	}
	@Override
	public T root(T t, int n) {
		throw new RuntimeException("Not implemented");
	}

        @Override
        public T sqrt(T t) {
          return root(t, 2);
        }
	
	@Override
	public T primitiveRoot() {
		if (!this.isFinite())
			throw new RuntimeException("Does not exist");
		if (this.primitiveRoot != null)
			return this.primitiveRoot;
		int q = this.getNumberOfElements() - 1;
		int t = q;
		List<Integer> factors = new ArrayList<Integer>();
		for (int f = 2; f <= t; f++) {
			if (t % f == 0) {
				factors.add(f);
				do {
					t /= f;
				} while(t % f == 0);
			}
		}
		T test;
		do {
			test = this.getRandomElement();
		} while (test.equals(this.zero()) || !this.isPrimitiveRoot(test, q, factors));
		this.primitiveRoot = test;
		return test;
	}
	private boolean isPrimitiveRoot(T t, int q, List<Integer> factors) {
		for (int f : factors) {
			if (this.power(t, q / f).equals(this.one()))
				return false;
		}
		return true;
	}
	@Override
	public boolean isIntegral() {
		return true;
	}
	@Override
	public boolean isEuclidean() {
		return true;
	}
	@Override
	public List<T> quotientAndRemainder(T dividend, T divisor) {
		List<T> result = new ArrayList<T>();
		result.add(this.divide(dividend, divisor));
		result.add(this.zero());
		return result;
	}

	@Override
	public int getNumberOfUnits() {
		return this.getNumberOfElements() - 1;
	}

	@Override
	public Iterable<T> getUnits() {
		return this.getNonZeroElements();
	}
	public T gcd(T t1, T t2) {
		return this.extendedEuclidean(t1, t2).get(0);
	}
	public List<T> extendedEuclidean(T t1, T t2) {
		List<T> result = new ArrayList<T>();
		if (t1.equals(this.zero()) && t2.equals(this.zero())) {
			result.add(this.zero());
			result.add(this.zero());
			result.add(this.zero());
			return result;
		}
		result.add(this.one());
		if (t1.equals(this.zero())) {
			result.add(this.zero());
			result.add(this.inverse(t2));
			return result;
		}
		if (t2.equals(this.zero()) || this.characteristic() == 2) {
			result.add(this.inverse(t1));
			result.add(this.zero());
			return result;
		}
		result.add(this.inverse(this.multiply(2, t1)));
		result.add(this.inverse(this.multiply(2, t2)));
		return result;
	}
}
