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

			@Override
			public boolean hasNext() {
				return true;
			}

			@Override
			public LexElement<T, S> next() {
				LexElement<T, S> next = i;
				increment(i);
				return next;
			}
		};
	}

	@Override
	public LexElement<T, S> least() {
		return new LexElement<>(first.least(), second.least());
	}
	
	@Override
	public LexElement<T, S> increment(LexElement<T, S> t) {
		return new LexElement<>(t.first, second.increment(t.second));
	}
}
