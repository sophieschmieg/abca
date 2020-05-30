package fields;

import java.math.BigInteger;
import java.util.Iterator;

public abstract class AbstractModule<T extends Element, S extends Element> implements Module<T, S> {

	@Override
	public Group<S> getAdditiveGroup() {
		return new Group<S>() {

			@Override
			public S getRandomElement() {
				return AbstractModule.this.getRandomElement();
			}

			@Override
			public boolean isFinite() {
				return AbstractModule.this.isFinite();
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				return AbstractModule.this.getNumberOfElements();
			}

			@Override
			public Iterator<S> iterator() {
				return AbstractModule.this.iterator();
			}

			@Override
			public S neutral() {
				return AbstractModule.this.zero();
			}

			@Override
			public S inverse(S s) {
				return AbstractModule.this.negative(s);
			}

			@Override
			public S operate(S s1, S s2) {
				return AbstractModule.this.add(s1, s2);
			}
			
		};
	}

	@Override
	public S add(S s1, S s2, S s3) {
		return add(s1, add(s2, s3));
	}

	@Override
	public S subtract(S minuend, S subtrahend) {
		return add(minuend, negative(subtrahend));
	}

	@Override
	public S scalarMultiply(int n, S s) {
		return scalarMultiply(getRing().getInteger(n), s);
	}

	@Override
	public S scalarMultiply(BigInteger n, S s) {
		return scalarMultiply(getRing().getInteger(n), s);
	}

	@Override
	public S scalarMultiply(T t1, T t2, S s) {
		return scalarMultiply(getRing().multiply(t1, t2), s);
	}

	@Override
	public S scalarMultiply(int n, T t, S s) {
		return scalarMultiply(getRing().multiply(n, t), s);
	}

	@Override
	public S scalarMultiply(BigInteger n, T t, S s) {
		return scalarMultiply(getRing().multiply(n, t), s);
	}

	@Override
	public S scalarMultiply(int n, T t1, T t2, S s) {
		return scalarMultiply(getRing().multiply(n, t1, t2), s);
	}

	@Override
	public S scalarMultiply(BigInteger n, T t1, T t2, S s) {
		return scalarMultiply(getRing().multiply(n, t1, t2), s);
	}

}
