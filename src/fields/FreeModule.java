package fields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class FreeModule<T extends Element> extends AbstractModule<T, Vector<T>> {
	private Ring<T> ring;
	private int dimension;

	public FreeModule(Ring<T> ring, int dimension) {
		this.ring = ring;
		this.dimension = dimension;
	}

	public Ring<T> getRing() {
		return ring;
	}

	@Override
	public Vector<T> zero() {
		List<T> coeffs = new ArrayList<T>();
		for (int i = 0; i < dimension; i++) {
			coeffs.add(ring.zero());
		}
		return new Vector<T>(coeffs);
	}

	@Override
	public Vector<T> add(Vector<T> s1, Vector<T> s2) {
		List<T> coeffs = new ArrayList<T>();
		for (int i = 1; i <= dimension; i++) {
			coeffs.add(ring.add(s1.get(i), s2.get(i)));
		}
		return new Vector<T>(coeffs);
	}

	@Override
	public Vector<T> negative(Vector<T> s) {
		List<T> coeffs = new ArrayList<T>();
		for (int i = 1; i <= dimension; i++) {
			coeffs.add(ring.negative(s.get(i)));
		}
		return new Vector<T>(coeffs);
	}

	@Override
	public Vector<T> scalarMultiply(T t, Vector<T> s) {
		List<T> coeffs = new ArrayList<T>();
		for (int i = 1; i <= dimension; i++) {
			coeffs.add(ring.multiply(t, s.get(i)));
		}
		return new Vector<T>(coeffs);
	}

	public T innerProduct(Vector<T> s1, Vector<T> s2) {
		T result = ring.zero();
		for (int i = 1; i <= dimension; i++) {
			result = ring.add(result, ring.multiply(s1.get(i), s2.get(i)));
		}
		return result;
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public Vector<T> getRandomElement() {
		List<T> coeffs = new ArrayList<T>();
		for (int i = 0; i < dimension; i++) {
			coeffs.add(ring.getRandomElement());
		}
		return new Vector<T>(coeffs);
	}

	@Override
	public boolean isFinite() {
		return ring.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return ring.getNumberOfElements().pow(dimension);
	}

	public List<Vector<T>> getGenerators() {
		List<Vector<T>> result = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			result.add(getUnitVector(i + 1));
		}
		return result;
	}

	public Vector<T> getUnitVector(int i) {
		List<T> coeffs = new ArrayList<>();
		for (int j = 1; j < i; j++) {
			coeffs.add(ring.zero());
		}
		coeffs.add(ring.one());
		for (int j = i + 1; j <= dimension; j++) {
			coeffs.add(ring.zero());
		}
		return new Vector<T>(coeffs);
	}

	@Override
	public Iterator<Vector<T>> iterator() {
		return new Iterator<Vector<T>>() {
			private List<Iterator<T>> it = null;
			private List<T> list = new ArrayList<T>();

			private void init() {
				if (this.it != null)
					return;
				this.it = new ArrayList<Iterator<T>>();
				for (int i = 0; i < dimension; i++) {
					this.it.add(ring.iterator());
					if (i == 0)
						this.list.add(null);
					else
						this.list.add(this.it.get(i).next());
				}
			}

			@Override
			public boolean hasNext() {
				init();
				for (int i = 0; i < dimension; i++) {
					if (this.it.get(i).hasNext())
						return true;
				}
				return false;
			}

			@Override
			public Vector<T> next() {
				init();
				boolean broken = false;
				for (int i = 0; i < dimension; i++) {
					if (this.it.get(i).hasNext()) {
						this.list.set(i, this.it.get(i).next());
						for (int j = 0; j < i; j++) {
							this.it.set(j, ring.iterator());
							this.list.set(j, this.it.get(j).next());
						}
						broken = true;
						break;
					}
				}
				if (!broken)
					throw new RuntimeException();
				return new Vector<T>(list);
			}
		};
	}
}
