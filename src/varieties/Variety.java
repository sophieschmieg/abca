package varieties;

import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.MathSet;

public interface Variety<T extends Element<T>> extends MathSet<ProjectivePoint<T>>{
	public Field<T> getField();
	public boolean hasRationalPoint(ProjectivePoint<T> p);
}
