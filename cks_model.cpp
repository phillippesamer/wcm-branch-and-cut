#include "cks_model.h"

CKSModel::CKSModel(IO *instance)
{

}

CKSModel::~CKSModel()
{

}

void CKSModel::create_variables()
{

}

void CKSModel::create_constraints()
{

}

void CKSModel::create_objective()
{

}

int CKSModel::solve(bool logging)
{
	return 0;
}

int CKSModel::save_optimization_status()
{
    return 0;
}

double CKSModel::runtime()
{
    return model->get(GRB_DoubleAttr_Runtime);
}

void CKSModel::set_time_limit(double tl)
{
    model->set(GRB_DoubleParam_TimeLimit, tl);
}

