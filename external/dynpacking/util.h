
const int data_pos = 60;

inline void date_time(ostream &o) {
	long t = time(0);
	o << ctime(&t) << endl;
}
inline void out_char (ostream &o, char ch, int count) {
	o << string (count, ch);
}
inline bool get_bool (istream &i, int pos = data_pos) {
	string line;
	getline (i, line);
	bool b;
	string data (line, pos);
	istringstream(data) >> b;
	return b;
}
inline int get_int (istream &i, int pos = data_pos) {
	string line;
	getline (i, line);
	int n;
	string data (line, pos);
	istringstream(data) >> n;
	return n;
}
inline double get_double (istream &i, int pos = data_pos) {
	string line;
	getline (i, line);
	double d;
	string data (line, pos);
	istringstream(data) >> d;
	return d;
}
inline void omit_line (istream &i) {
	string line;
	getline (i, line);
}
inline void ivar(ostream &o, string sname_p, int kvalue_p, string smean_p) {
	o.setf(ios::left, ios::adjustfield);
	o << setw (12) << sname_p;
	o.setf(ios::right, ios::adjustfield);
	o << setw (13) << kvalue_p;
	o.fill ('.');
	o << setw (45) << smean_p;
	o.fill (' ');
	o << endl;
}
inline void lvar(ostream &o, string sname_p, long kvalue_p, string smean_p) {
	o.setf(ios::left, ios::adjustfield);
	o << setw (12) << sname_p;
	o.setf(ios::right, ios::adjustfield);
	o << setw (13) << kvalue_p;
	o.fill ('.');
	o << setw (45) << smean_p;
	o.fill (' ');
	o << endl;
}
inline void rvar(ostream &o, string sname_p, double dvalue_p, string smean_p) {
	o.setf(ios::left, ios::adjustfield);
	o << setw (12) << sname_p;
	o.setf(ios::right, ios::adjustfield);
	o << setw (13) << dvalue_p;
	o.fill ('.');
	o << setw (45) << smean_p;
	o.fill (' ');
	o << endl;
}
inline void open_failure (ostream &o, const char *s) {
	o << "Cannot open " << s << " file" << endl;
	exit(1);
}
inline void error (ostream &o, const char *s) {
	o << "Error: " << s << endl;
	exit(1);
}
inline void blines (ostream &o, int n = 1) {
	for (int i = 0; i < n; i++) o << endl;
}
inline void page (ostream &o) {
	o << endl << (char) 12 << endl;	
}

inline double max (double a, double b) {
	return (a >= b) ? a : b;
}

inline double min (double a, double b) {
	return (a <= b) ? a : b;
}

inline double sign (double a, double b) {
	return (b >= 0) ? fabs(a) : -fabs(a);
}

