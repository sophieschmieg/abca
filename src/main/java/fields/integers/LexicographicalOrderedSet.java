package fields.integers;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.integers.LexicographicalOrderedSet.LexElement;
import fields.interfaces.Element;
import fields.interfaces.WellOrder;

public class LexicographicalOrderedSet<T extends Element<T>, S extends Element<S>>
		implements WellOrder<LexElement<T, S>> {
	private WellOrder<T> first;
	private WellOrder<S> second;

	public static class LexElement<T extends Element<T>, S extends Element<S>>
			extends AbstractElement<LexElement<T, S>> {
		private T first;
		private S second;

		public LexElement(T first, S second) {
			this.first = first;
			this.second = second;
		}

		@Override
		public int compareTo(LexElement<T, S> o) {
			int cmp = first.compareTo(o.first);
			if (cmp != 0) {
				return cmp;
			}
			return second.compareTo(o.second);
		}

		public String toString() {
			return "(" + first.toString() + ", " + second.toString() + ")";
		}
	}

	public LexicographicalOrderedSet(WellOrder<T> first, WellOrder<S> second) {
		this.first = first;
		this.second = second;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public LexElement<T, S> getRandomElement() {
		return new LexElement<>(first.getRandomElement(), second.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return first.isFinite() && second.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return first.getNumberOfElements().multiply(second.getNumberOfElements());
	}

	@Override
	public Iterator<LexElement<T, S>> iterator() {
		return new Iterator<>() {
			LexElement<T, S> i = least();
			Iterator<T> firstIt = first.iterator();
			Iterator<S> secondIt = second.iterator();

			@Override
			public boolean hasNext() {
				return firstIt.hasNext();
			}

			@Override
			public LexElement<T, S> next() {
				LexElement<T, S> next = i;
				if (secondIt.hasNext()) {
					i = new LexElement<>(next.first, secondIt.next());
				} else {
					secondIt = second.iterator();
					i = new LexElement<>(firstIt.next(), secondIt.next());
				}
				return next;
			}
		};
	}

	@Override
	public LexElement<T, S> least() {
		return new LexElement<>(first.least(), second.least());
	}

	@Override
	public int compare(LexElement<T, S> t1, LexElement<T, S> t2) {
		int cmp = t1.first.compareTo(t2.first);
		if (cmp != 0) {
			return cmp;
		}
		return t1.second.compareTo(t2.second);
	}
}
