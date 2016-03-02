
// >>>>>>>>>> DEFINITIONS <<<<<<<<<<
#define slices 50
#define stacks 190

class cRGB_Byte_Pixel
{
public:
	cRGB_Byte_Pixel()
	{
		r=255;
		g=255;
		b=255;
		a=255;
	};

	union
	{
		struct
		{
			int r;
			int g;
			int b;
			int a;
		};
		int rgb[4];
	};
}; 

GLuint Tex1, Tex2, Tex3, Tex4, Tex5, Tex6, Tex7, Tex8, Tex9, Tex10, Tex11, Tex12, Tex13, Tex14, Tex15, TexCue;
GLuint List1, List2, List3, List4, List5, List6, List7, List8, List9, List10, List11, List12, List13, List14, List15, ListCue;
int TexWidth, TexHeight;
cRGB_Byte_Pixel * TD1, * TD2, * TD3, * TD4, * TD5, * TD6, * TD7, * TD8, * TD9, * TD10, * TD11, * TD12, * TD13, * TD14, * TD15, * TDCue;

// >>>>>>>>>> FUNCTIONS <<<<<<<<<<
void CreateSphere(double, int, int, GLuint);
void DrawSphere(GLuint);
