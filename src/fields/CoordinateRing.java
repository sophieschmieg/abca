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

public class CoordinateRing<T extends Element> extends AbstractRing<CoordinateRingElement<T>> implements Ring<CoordinateRingElement<T>> {
	private PolynomialRing<T> ring;
        private PolynomialRing<T>.Ideal ideal;
	
	public CoordinateRing(PolynomialRing<T> ring, PolynomialRing<T>.Ideal ideal) {
          this.ring = ring;
          this.ideal = ideal;
	}
	@Override
	public boolean isFinite() {
		return this.ring.getField().isFinite() && this.ideal.dimension() == 0;
	}
        @Override
        public String toString() {
          return this.ring.toString() + "/" + this.ideal.toString();
        }
	public CoordinateRingElement<T> getEmbedding(Polynomial<T> t) {
		return new CoordinateRingElement(this.ring, this.ideal, t);
	}
	@Override
	public ExtensionFieldElement<T> zero() {
		return this.getEmbedding(this.ring.zero());
	}
	@Override
	public ExtensionFieldElement<T> one() {
		return this.getEmbedding(this.ring.one());
	}
	@Override
	public int characteristic() {
		return this.ring.characteristic();
	}
	@Override
	public CoordinateRingElement<T> add(CoordinateRingElement<T> t1,
			CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.add(t1.polynomial, t2.polynomial));
	}
	@Override
	public CoordinateRingElement<T> negative(CoordinateRingElement<T> t) {
		return this.getEmbedding(this.ring.negative(t.polynomial));
	}
	@Override
	public CoordinateRingElement<T> multiply(CoordinateRingElement<T> t1,
			CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1.polynomial, t2.polynomial));
	}
	@Override
	public CoordinateRingElement<T> inverse(CoordinateRingElement<T> t) {
		return this.getEmbedding(this.ring.inverse(t.polynomial));
	}
        @Override
        public boolean isUnit(CoordinateRingElement<T> t) {
          return this.ring.isUnit(t.polynomial);
        }
	@Override
	public ExtensionFieldElement<T> getRandomElement() {
		throw new UnsupportedOperationException();
	}
	@Override
	public int getNumberOfElements() throws InfinityException {
		return -1;
	}
	@Override
	public Iterable<ExtensionFieldElement<T>> getElements()
			throws InfinityException {
		throw new UnsupportedOperationException();
	}
	public static class CoordinateRingElement<T extends Element> implements Element  {
		private PolynomialRing<T> ring;
		private PolynomialRing<T>.Ideal ideal;
                private Polynomial<T> polynomial;
		
		public CoordinateRingElement(PolynomialRing<T> ring, PolynomialRing<T>.Ideal ideal, Polynomial<T> polynomial) {
	          this.ring = ring;
                  this.ideal = ideal;
                  this.polynomial = ideal.residue(polynomial);
                }
		public boolean equals(Object o) {
			if (!(o instanceof CoordinateRing.CoordinateRingElement))
				return false;
			@SuppressWarnings("unchecked")
			CoordinateRingElement<T> t = (CoordinateRingElement<T>)o;
			return this.polynomial.equals(t.polynomial);
		}
		public String toString() {
                  return this.polynomial.toString();
		}
		public List<T> getElement() {
			return Collections.unmodifiableList(this.coeff);
		}
		@Override
		public int compareTo(Element e) {
			@SuppressWarnings("unchecked")
			CoordinateRingElement<T> o = (CoordinateRingElement<T>)e;
			return this.polynomial.compareTo(o.polynomial);
		}
	}
}
