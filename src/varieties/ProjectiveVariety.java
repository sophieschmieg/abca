package varieties;

import fields.Element;
import fields.Field;
import fields.InfinityException;
import fields.PolynomialRing;

public class ProjectiveVariety<T extends Element> implements Variety<T> {
	private ProjectiveSpace<T> space;
	private PolynomialRing<T>.Ideal definingIdeal; 

	public ProjectiveVariety(ProjectiveSpace<T> space, PolynomialRing<T>.Ideal definingIdeal) {
		this.space = space;
		this.definingIdeal = definingIdeal;
	}
	
	public ProjectiveSpace<T> getSpace() {
		return this.space;
	}
	@Override
	public ProjectivePoint<T> getRandomElement() {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isFinite() {
		return this.space.isFinite() || this.definingIdeal.dimension() == 0;
	}

	@Override
	public int getNumberOfElements() throws InfinityException {
		throw new UnsupportedOperationException();
	}

	@Override
	public Iterable<ProjectivePoint<T>> getElements() throws InfinityException {
		throw new UnsupportedOperationException();
	}

	@Override
	public Field<T> getField() {
		return this.space.getField();
	}

	@Override
	public boolean hasRationalPoint(ProjectivePoint<T> p) {
		return this.space.asIdeal(p).contains(this.definingIdeal);
	}
	public boolean contains(ProjectiveVariety<T> var) {
		return var.definingIdeal.contains(this.definingIdeal);
	}
	public ProjectiveVariety<T> intersect(ProjectiveVariety<T> var) {
		return new ProjectiveVariety<T>(this.space, this.definingIdeal.add(var.definingIdeal));
	}
	public ProjectiveVariety<T> unite(ProjectiveVariety<T> var) {
		return new ProjectiveVariety<T>(this.space, this.definingIdeal.intersect(var.definingIdeal));
	}
}
