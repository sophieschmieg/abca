package fields.helper;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Ring;
import fields.local.Value;
import fields.vectors.Vector;

public abstract class AbstractIdeal<T extends Element<T>> extends AbstractModule<T, T> implements Ideal<T> {
	private Ring<T> r;

	public AbstractIdeal(Ring<T> r) {
		this.r = r;
	}

	@Override
	public Exactness exactness() {
		return r.exactness();
	}

	@Override
	public Ring<T> getRing() {
		return r;
	}
	
	@Override
	public boolean isRadical() {
		return r.radical(this).equals(this);
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
		return r.multiply(s, t);
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
		for (T g : this.generators()) {
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
				for (int i = 0; i < size; i++) {
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
		if (s.isEmpty()) {
			return true;
		}
		List<T> generators = generators();
		if (generators.isEmpty() ||( generators.size() == 1 && generators.get(0).equals(r.zero()))) {
			return true;
		}
		if (s.size() > 1) {
			return false;
		}
		return r.isZeroDivisor(s.get(0));
	}

	@Override
	public final List<List<T>> nonTrivialCombinations(List<T> s) {
		if (s.size() > 1) {
			List<List<T>> result = new ArrayList<>();
			for (int i = 1; i < s.size(); i++) {
				List<T> here = new ArrayList<>();
				here.add(s.get(i));
				for (int j = 2; j < i; j++) {
					here.add(r.zero());
				}
				here.add(r.negative(s.get(0)));
				for (int j = i + 1; j < s.size(); j++) {
					here.add(r.zero());
				}
				result.add(here);
			}
			return result;
		}
		if (s.isEmpty() || !r.isZeroDivisor(s.get(0))) {
			throw new ArithmeticException("List is independent");
		}
		if (s.get(0).equals(r.zero())) {
			return Collections.singletonList(Collections.singletonList(r.one()));
		}
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean contains(Ideal<T> other) {
		for (T t : other.generators()) {
			if (!contains(t)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean equalsIdeal(Ideal<T> other) {
		return contains(other) && other.contains(this);
	}

	@SuppressWarnings("unchecked")
	@Override
	public final boolean equals(Object other) {
		return equalsIdeal((Ideal<T>) other);
	}

	@Override
	public int compareTo(Ideal<T> other) {
		if (equalsIdeal(other)) {
			return 0;
		}
		if (contains(other)) {
			return 1;
		}
		if (other.contains(this)) {
			return -1;
		}
		List<T> generatorsThis = new ArrayList<>();
		generatorsThis.addAll(generators());
		List<T> generatorsOther = new ArrayList<>();
		generatorsOther.addAll(other.generators());
		if (generatorsThis.size() != generatorsOther.size()) {
			return generatorsThis.size() - generatorsOther.size();
		}
		generatorsThis.sort(null);
		generatorsOther.sort(null);
		for (int i = 0; i < generatorsThis.size(); i++) {
			int cmp = generatorsThis.get(i).compareTo(generatorsOther.get(i));
			if (cmp != 0) {
				return 0;
			}
		}
		throw new RuntimeException("equalsIdeal was wrong");
	}

	@Override
	public final Vector<T> asVector(T s) {
		return new Vector<>(generate(s));
	}

	@Override
	public final boolean isGeneratingModule(List<T> s) {
		return equalsIdeal(r.getIdeal(s));
	}

	@Override
	public Value maximumPowerContains(T t) {
		if (t.equals(r.zero()) || contains(r.one())) {
			return Value.INFINITY;
		}
		int value = 0;
		Ideal<T> ideal = this;
		while (ideal.contains(t)) {
			value++;
			ideal = r.multiply(ideal, this);
		}
		return new Value(value);
	}

	@Override
	public String toString() {
		StringBuilder build = new StringBuilder();
		build.append("(");
		boolean first = true;
		for (T g : generators()) {
			if (first) {
				first = false;
			} else {
				build.append(", ");
			}
			build.append(g.toString());
		}
		build.append(")");
		return build.toString();
	}

	@Override
	public boolean isLeftIdeal() {
		if (r.isCommutative()) {
			return true;
		}
		throw new UnsupportedOperationException("Not implemented");
	}

	@Override
	public boolean isRightIdeal() {
		if (r.isCommutative()) {
			return true;
		}
		throw new UnsupportedOperationException("Not implemented");
	}
}
