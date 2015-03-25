#ifndef PCA_H
#define PCA_H

#include <iostream>
#include <vector>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

class PCA
{
public:
  static std::vector<double> compute( Eigen::MatrixXd D )
  {
    std::vector<double> eigenvalues;
    std::vector< std::vector<double> > eigenvectors;
    
    std::size_t m = D.rows();
    std::size_t n = D.cols();
    double mean;
    Eigen::VectorXd meanVector;
    
    for(std::size_t i = 0; i < n; i++)
    {
      // compute mean
      mean = (D.col(i).sum())/(double)m;
      // create a vector with constant value = mean
      meanVector = Eigen::VectorXd::Constant(m, mean); 
      D.col(i) -= meanVector;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd( D, Eigen::ComputeThinV );
    
    // eigenvalues
    Eigen::JacobiSVD<Eigen::MatrixXd>::SingularValuesType singularValues = svd.singularValues();  
    Eigen::VectorXd ev = Eigen::VectorXd::Zero( m );
    for( std::size_t i = 0; i < singularValues.size(); i++ )
    {
      double s = static_cast<double>( singularValues(i) );
      eigenvalues.push_back( s * s / static_cast<double>( m - 1 ) );
    }
    
    // eigenvectors
    std::vector< std::vector<double> >().swap( eigenvectors );
    eigenvectors.resize( n, std::vector<double>( n ) );
    Eigen::MatrixXd mat = svd.matrixV();
    for( std::size_t index = 0; index < n; index++ )
      for( std::size_t d = 0; d < n; d++ )
        eigenvectors[index][d] = mat(d,index);
    
    /*
    for( std::size_t i = 0; i < eigenvalues.size(); i++ )
      std::cout << eigenvalues.at(i) << std::endl;
    
    for( std::size_t i = 0; i < eigenvectors.size(); i++ )
    {
      for( std::size_t j = 0; j < eigenvectors.at(i).size(); j++ )
        std::cout << eigenvectors.at(i).at(j) << std::endl;
      
      std::cout << std::endl;
    }
    */
    
    std::vector<double> longestPrincipalComponent;
    longestPrincipalComponent.resize(3);
    for( std::size_t i = 0; i < eigenvectors.front().size(); i++ )
      longestPrincipalComponent.at(i) = eigenvectors.at(i).at(0);
    
    return longestPrincipalComponent;
  }
};

#endif // PCA_H