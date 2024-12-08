#ifndef MOLECULE_DATA_H
#define MOLECULE_DATA_H

typedef enum
{
    LI_2_SINGLET, LI_2_TRIPLET,
    NA_2_SINGLET, NA_2_TRIPLET,
    K_2_SINGLET, K_2_TRIPLET,
    RB_2_SINGLET, RB_2_TRIPLET
} target;


double get_polarizabiliy_anisotropy(target species, )
{
    switch(species)
    {
    case LI_2_SINGLET:
        return 0.0;
        break;

    case LI_2_TRIPLET:
        return ;
        break;

    case NA_2_SINGLET: 
        	return
    }
}


#endif 