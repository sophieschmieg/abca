package varieties;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.helper.CoordinateRing;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;

public class ProjectiveVariety<T extends Element<T>> implements Variety<T> {
	private ProjectiveSpace<T> space;
	private Ideal<Polynomial<T>> definingIdeal;
	private CoordinateRing<T> coordinateRing;

	public ProjectiveVariety(ProjectiveSpace<T> space, Ideal<Polynomial<T>> definingIdeal) {
		this.space = space;
		this.definingIdeal = definingIdeal;
		this.coordinateRing = new CoordinateRing<T>(space.getRing(), this.definingIdeal);
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
		return this.space.isFinite() || this.coordinateRing.krullDimension() == 0;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new UnsupportedOperationException();
	}

	@Override
	public Iterator<ProjectivePoint<T>> iterator() {
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
		return new ProjectiveVariety<T>(this.space, this.space.getRing().add(this.definingIdeal, var.definingIdeal));
	}

	public ProjectiveVariety<T> unite(ProjectiveVariety<T> var) {
		return new ProjectiveVariety<T>(this.space,
				this.space.getRing().intersect(this.definingIdeal, var.definingIdeal));
	}
}
