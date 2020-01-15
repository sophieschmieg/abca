package fields;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.ExtensionField.ExtensionFieldElement;
import fields.Polynomial.Monomial;

public class ExtensionField<T extends Element> extends AbstractField<ExtensionFieldElement<T>> implements Field<ExtensionFieldElement<T>> {
	private Polynomial<T> minimalpolynomial;
	private Field<T> field;
	private int degree;
	private PolynomialRing<T> ring;
	private Map<ExtensionFieldElement<T>, ExtensionFieldElement<T>> inverses;
	private Map<T, ExtensionFieldElement<T>> embed;
	private Comparator<Polynomial.Monomial> comparator;
	
	public ExtensionField(Polynomial<T> minimalpolynomial, Field<T> field) {
		this.comparator = new Polynomial.LexographicalOrder();
		this.field = field;
		this.ring = minimalpolynomial.getRing();
		if (minimalpolynomial.getNumVars() != 1) // TODO: Irreducible
			throw new ArithmeticException("Not a minimal polynomial");
		this.degree = minimalpolynomial.getDegree();
		this.minimalpolynomial = minimalpolynomial.normalize();
		this.inverses = new TreeMap<ExtensionFieldElement<T>, ExtensionFieldElement<T>>();
		this.embed = new TreeMap<T, ExtensionFieldElement<T>>();
	}
	@Override
	public boolean isFinite() {
		return this.field.isFinite();
	}
        @Override
        public String toString() {
          return this.field.toString() + "[X]/(" + this.minimalpolynomial.toString() + ")";
        }
	public ExtensionFieldElement<T> getEmbedding(T t) {
		if (this.embed.containsKey(t))
			return this.embed.get(t);
		List<T> list = new ArrayList<T>();
		list.add(t);
		for (int i = 1; i < this.degree; i++)
			list.add(this.field.zero());
		ExtensionFieldElement<T> e = new ExtensionFieldElement<T>(this.field, this.minimalpolynomial, list);
		this.embed.put(t, e);
		return this.getEmbedding(t);
	}
	@Override
	public ExtensionFieldElement<T> zero() {
		return this.getEmbedding(this.field.zero());
	}
	@Override
	public ExtensionFieldElement<T> one() {
		return this.getEmbedding(this.field.one());
	}
	@Override
	public int characteristic() {
		return this.field.characteristic();
	}
	@Override
	public ExtensionFieldElement<T> add(ExtensionFieldElement<T> t1,
			ExtensionFieldElement<T> t2) {
		List<T> list = new ArrayList<T>();
		for (int i = 0; i < this.degree; i++)
			list.add(this.field.add(t1.getElement().get(i), t2.getElement().get(i)));
		return new ExtensionFieldElement<T>(this.field, this.minimalpolynomial, list);
	}
	@Override
	public ExtensionFieldElement<T> negative(ExtensionFieldElement<T> t) {
		List<T> list = new ArrayList<T>();
		for (int i = 0; i < this.degree; i++)
			list.add(this.field.negative(t.getElement().get(i)));
		return new ExtensionFieldElement<T>(this.field, this.minimalpolynomial, list);
	}
	@Override
	public ExtensionFieldElement<T> multiply(ExtensionFieldElement<T> t1,
			ExtensionFieldElement<T> t2) {
		Polynomial<T> res = this.ring.multiply(this.asPolynomial(t1), this.asPolynomial(t2));
		return this.fromPolynomial(this.ring.quotientAndRemainder(res, this.minimalpolynomial).get(1));
	}
	@Override
	public ExtensionFieldElement<T> inverse(ExtensionFieldElement<T> t) {
          List<Polynomial<T>> egcd = this.ring.extendedEuclidean(this.asPolynomial(t), this.minimalpolynomial);
          if (!this.ring.isUnit(egcd.get(0))) {
            throw new ArithmeticException("Division by zero");
          }
          ExtensionFieldElement<T> inv = this.fromPolynomial(this.ring.multiply(egcd.get(1), this.ring.inverse(egcd.get(0))));
		return inv;
	}
	private Polynomial<T> asPolynomial(ExtensionFieldElement<T> t) {
		Map<Monomial,T> map = new TreeMap<Monomial, T>(this.comparator);
		for (int i = 0; i < this.degree; i++)
			map.put(new Monomial(new int[] {i}), t.getElement().get(i));
		return new Polynomial<T>(map, this.ring);
	}
	private ExtensionFieldElement<T> fromPolynomial(Polynomial<T> t) {
		List<T> list = new ArrayList<T>();
		for (int i = 0; i < this.degree; i++)
			list.add(t.getCoefficient(new Monomial(new int[] {i})));
		return new ExtensionFieldElement<T>(this.field, this.minimalpolynomial, list);
	}
	@Override
	public ExtensionFieldElement<T> getRandomElement() {
		List<T> list = new ArrayList<T>();
		for (int i = 0; i < this.degree; i++)
			list.add(this.field.getRandomElement());
		return new ExtensionFieldElement<T>(this.field, this.minimalpolynomial, list);
	}
	@Override
	public int getNumberOfElements() throws InfinityException {
		int num = this.field.getNumberOfElements();
		int res = 1;
		for (int i = 0; i < this.degree; i++)
			res *= num;
		
		return res;
	}
	@Override
	public Iterable<ExtensionFieldElement<T>> getElements()
			throws InfinityException {
		final Iterable<T> elements = this.field.getElements();
		return new Iterable<ExtensionFieldElement<T>>() {
			@Override
			public Iterator<ExtensionFieldElement<T>> iterator() {
				return new Iterator<ExtensionFieldElement<T>>() {
					private List<Iterator<T>> it = null;
					private List<T> list = new ArrayList<T>();
					
					private void init() {
						if (this.it != null)
							return;
						this.it = new ArrayList<Iterator<T>>();
						for (int i = 0; i < degree; i++) {
							this.it.add(elements.iterator());
							if (i == 0)
								this.list.add(null);
							else
								this.list.add(this.it.get(i).next());
						}
					}
					
					@Override
					public boolean hasNext() {
						init();
						for (int i = 0; i < degree; i++) {
							if (this.it.get(i).hasNext())
								return true;
						}
						return false;
					}

					@Override
					public ExtensionFieldElement<T> next() {
						init();
						boolean broken = false;
						for (int i = 0; i < degree; i++) {
							if (this.it.get(i).hasNext()) {
								this.list.remove(i);
								this.list.add(i, this.it.get(i).next());
								for (int j = 0; j < i; j++) {
									this.it.remove(j);
									this.it.add(j, elements.iterator());
									this.list.remove(j);
									this.list.add(j, this.it.get(j).next());
								}
								broken = true;
								break;
							}
						}
						if (!broken)
							throw new RuntimeException();
						return new ExtensionFieldElement<T>(field, minimalpolynomial, list);
					}

					@Override
					public void remove() {
						throw new UnsupportedOperationException();
					}
				};
			}
		};
	}
	public static class ExtensionFieldElement<T extends Element> implements Element  {
		private Field<T> base;
		private List<T> coeff;
		
		public ExtensionFieldElement(Field<T> base, Polynomial<T> minimalpolynomial, List<T> coeff) {
			this.base = base;
			this.coeff = new ArrayList<T>();
			if (minimalpolynomial.getNumVars() != 1 || coeff.size() != minimalpolynomial.getDegree())
				throw new ArithmeticException("Invalid minimal polynomial");
			this.coeff.addAll(coeff);
		}
		public boolean equals(Object o) {
			if (!(o instanceof ExtensionField.ExtensionFieldElement))
				return false;
			@SuppressWarnings("unchecked")
			ExtensionFieldElement<T> t = (ExtensionFieldElement<T>)o;
			return this.coeff.equals(t.coeff);
		}
		public String toString() {
			boolean first = true;
			boolean onlyzero = true;
			StringBuffer buf = new StringBuffer();
			for (int i = 0; i < this.coeff.size(); i++) {
				if (!this.coeff.get(i).equals(this.base.zero())) {
					onlyzero = false;
					if (first)
						first = false;
					else
						buf.append(" + ");
					if (!this.coeff.get(i).equals(this.base.one()) || i == 0)
						buf.append(this.coeff.get(i));
					if (i == 1)
						buf.append("a");
					if (i > 1)
						buf.append("a^" + i);
				}
			}
			if (onlyzero)
				return "0";
			return buf.toString();
		}
		public List<T> getElement() {
			return Collections.unmodifiableList(this.coeff);
		}
		@Override
		public int compareTo(Element e) {
			@SuppressWarnings("unchecked")
			ExtensionFieldElement<T> o = (ExtensionFieldElement<T>)e;
			for (int i = 0; i < this.coeff.size(); i++) {
				int cmp = this.coeff.get(i).compareTo(o.coeff.get(i));
				if (cmp != 0)
					return cmp;
			}
			return 0;
		}
	}
	public static<T extends Element> PolynomialRing<T> getPolynomialRing(Field<T> base) {
		return new PolynomialRing<T>(base, 1, Polynomial.LEX);
	}
}
