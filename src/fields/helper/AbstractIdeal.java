package fields.helper;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Ring;
import fields.vectors.Vector;

public abstract class AbstractIdeal<T extends Element<T>> extends AbstractModule<T, T> implements Ideal<T> {
	private Ring<T> r;

	public AbstractIdeal(Ring<T> r) {
		this.r = r;
	}

	@Override
	public Ring<T> getRing() {
		return r;
	}

	@Override
	public T zero() {
		return r.zero();
	}

	@Override
	public T add(T s1, T s2) {
		return r.add(s1, s2);
	}

	@Override
	public T negative(T s) {
		return r.negative(s);
	}

	@Override
	public T scalarMultiply(T t, T s) {
		return r.multiply(t, s);
	}

	public final boolean isPrincipal() {
		return this.generators().size() == 1;
	}
	
	@Override
	public final boolean isFree() {
		return isPrincipal();
	}

	@Override
	public final T getRandomElement() {
		T result = r.zero();
		for (T g :this.generators()) {
			result = r.add(result, r.multiply(r.getRandomElement(), g));
		}
		return result;
	}

	@Override
	public final Iterator<T> iterator() {
		return new Iterator<T>() {
			private List<Iterator<T>> it = null;
			private List<T> list = new ArrayList<T>();
			private int size = generators().size();

			private void init() {
				if (this.it != null)
					return;
				this.it = new ArrayList<Iterator<T>>();
				for (int i = 0; i < size; i++) {
					this.it.add(r.iterator());
					if (i == 0)
						this.list.add(null);
					else
						this.list.add(this.it.get(i).next());
				}
			}

			@Override
			public boolean hasNext() {
				init();
				for (int i = 0; i < size; i++) {
					if (this.it.get(i).hasNext())
						return true;
				}
				return false;
			}

			@Override
			public T next() {
				init();
				boolean broken = false;
				for (int i = 0; i < size; i++) {
					if (this.it.get(i).hasNext()) {
						this.list.set(i, this.it.get(i).next());
						for (int j = 0; j < i; j++) {
							this.it.set(j, r.iterator());
							this.list.set(j, this.it.get(j).next());
						}
						broken = true;
						break;
					}
				}
				if (!broken)
					throw new RuntimeException();
				T result = r.zero();
				List<T> g = generators();
				for (int i = 0; i<size;i++) {
					result = r.add(result, r.multiply(list.get(i), g.get(i)));
				}
				return result;
			}
		};
	}
	
	@Override
	public final List<T> getModuleGenerators() {
		return generators();
	}
	
	@Override
	public final boolean isLinearIndependent(List<T> s) {
		return s.size() == 0;
	}
	
	@Override
	public final boolean contains(Ideal<T> other) {
		for (T t : other.generators()) {
			if (!contains(t)) {
				return false;
			}
		}
		return true;
	}
	
	@Override
	public final boolean equalsIdeal(Ideal<T> other) {
		return contains(other) && other.contains(this);
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public final boolean equals(Object other) {
		return equalsIdeal((Ideal<T>)other);
	}
	
	@Override
	public final Vector<T> asVector(T s) {
		return new Vector<>(generate(s));
	}
	
	@Override
	public final boolean isGeneratingModule(List<T> s) {
		return equalsIdeal(r.getIdeal(s));
	}
}
