#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
//#include "itkAffineTransform.h"
#include "itkTranslationTransform.h"
#include <itkMatrix.h>
#include <itkExtractImageFilter.h>

#include "igtlOSUtil.h"
#include "igtlMessageHeader.h"
#include "igtlTransformMessage.h"
#include "igtlPositionMessage.h"
#include "igtlImageMessage.h"
#include "igtlClientSocket.h"
#include "igtlStatusMessage.h"
#include <igtl_util.h>

const    unsigned int    DimensionImage = 2;
const    unsigned int    DimensionVolume = 3;
typedef  unsigned char   PixelType;
typedef  itk::Image< PixelType, DimensionImage >  ImageType;
typedef  itk::Image< PixelType, DimensionVolume >  VolumeType;
//typedef itk::AffineTransform< double, DimensionImage > TransformType;
typedef itk::TranslationTransform< double, DimensionImage > TransformType;
typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
typedef itk::MeanSquaresImageToImageMetric< ImageType, ImageType > MetricType;
typedef itk::LinearInterpolateImageFunction< ImageType, double > InterpolatorType;
typedef itk::ImageRegistrationMethod< ImageType, ImageType > RegistrationType;
typedef itk::ImageFileReader<VolumeType> FileReaderType;
typedef itk::ExtractImageFilter< VolumeType, ImageType > FilterType;
typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
typedef RegistrationType::ParametersType ParametersType;


class Images
{
  public:

  Images();
  int IGTtoITKImage();
  int ITKtoIGTImage();
  void SetParametersFromIGT();
  void SetParametersFromITK();

  ImageType::Pointer imageData;
  igtl::ImageMessage::Pointer imgMsg;

  protected:

  int sizeIm[3];
  float originIm[3];
  float spacingIm[3];

};

Images::Images()
{

  // Create imageMessage and ITKimage for this image
  imageData  = ImageType::New();
  imgMsg = igtl::ImageMessage::New();

}

void Images::SetParametersFromIGT()
{

  this->imgMsg->GetDimensions(this->sizeIm);
  this->imgMsg->GetOrigin(this->originIm);
  this->imgMsg->GetSpacing(this->spacingIm);

  return;

}

void Images::SetParametersFromITK()
{

  ImageType::SizeType size = this->imageData->GetLargestPossibleRegion().GetSize();
  sizeIm[0] = size[0]; sizeIm[1] = size[1]; sizeIm[2] = 1;
  ImageType::PointType origin = this->imageData->GetOrigin();
  originIm[0] = origin[0]; originIm[1] = origin[1]; originIm[2] = origin[2];
  ImageType::SpacingType spacing = this->imageData->GetSpacing();
  spacingIm[0] = spacing[0]; spacingIm[1] = spacing[1]; spacingIm[2] = spacing[2];

  return;

}

int Images::IGTtoITKImage()
{

  // Retrieve the image data from image message
  int       size[3];          // image dimension (pixels)
  float     spacing[3];       // spacing (mm/pixel)
  float     spacingITK[0];    // spacing for 2D image (mm/pixel)
  int       endian;           // endian (not used)
  float     origin[3];        // origin ()
  igtl::Matrix4x4   matrix;   // image origin and orientation matrix

  endian = this->imgMsg->GetEndian();
  this->imgMsg->GetDimensions( size );
  this->imgMsg->GetSpacing( spacing );
  this->imgMsg->GetMatrix( matrix );
  this->imgMsg->GetOrigin( origin );

  // Set image data to ITK image (3D information -> 2D)
  ImageType::RegionType     region;
  ImageType::IndexType      start;
  start[0] = origin[0];     start[1] = origin[1];
  ImageType::SizeType       sizeregion;
  sizeregion[0] = size[0];  sizeregion[1] = size[1];
  region.SetSize( sizeregion );
  region.SetIndex( start );
  spacingITK[0] = spacing[0]; spacingITK[1] = spacing[1];

  this->imageData->SetRegions( region );
  this->imageData->SetSpacing( spacingITK );
  this->imageData->Allocate();

  // Copy image data into ITK image
  memcpy( this->imageData->GetBufferPointer(), this->imgMsg->GetScalarPointer(), this->imgMsg->GetSubVolumeImageSize() );
  this->imageData->Modified();

  return 1;

}

int Images::ITKtoIGTImage()
{

  // Retrieve the image data from ITK image
  ImageType::SizeType       size;
  ImageType::IndexType      start;
  ImageType::SpacingType    spacing;

  size = this->imageData->GetLargestPossibleRegion().GetSize();
  start = this->imageData->GetLargestPossibleRegion().GetIndex();
  spacing = this->imageData->GetSpacing();

  // Set image data to image message (2D information -> 3D)
  int                       sizeIGT[3];
  ImageType::RegionType     region;
  float                     spacingIGT[3];  // spacing (mm/pixel)
  int                       scalarType;     // always UINT8
  //int                     svoffset[3];    // sub-volume offset
  //int                     svsize[3];      // sub-volume size
  sizeIGT[0] = size[0];         sizeIGT[1] = size[1];           sizeIGT[2] = 1;
  spacingIGT[0] = spacing[0];   spacingIGT[1] = spacing[1];     spacingIGT[2] = spacing[1]; // third dimension?
  //svsize[0] = sizeIGT[0];     svsize[1] = sizeIGT[1];         svsize[2] = sizeIGT[2];
  //svoffset[0] = start[0];     svoffset[1] = start[1];         svoffset[2] = start[2];
  scalarType = igtl::ImageMessage::TYPE_UINT8;

  this->imgMsg->SetDimensions( sizeIGT );
  this->imgMsg->SetSpacing( spacingIGT );
  this->imgMsg->SetScalarType( scalarType );
  //this->imgMsg->SetSubVolume( sizeIGT, svoffset );
  this->imgMsg->AllocateScalars();

  // Copy image data into ITK image
  memcpy( this->imgMsg->GetScalarPointer(), imageData->GetBufferPointer(), this->imgMsg->GetSubVolumeImageSize() );
  // Pack (serialize) and send
  this->imgMsg->Pack();

  return 1;

}

class Volumes
{

  public:

  Volumes();
  int ITKtoIGTVolume();

  VolumeType::Pointer volumeData;
  igtl::ImageMessage::Pointer imgMsg;

  protected:

  int sizeVol[3];
  float originVol[3];
  float spacingVol[3];

};

Volumes::Volumes()
{

  // Create imageMessage and ITKvolume for this volume
  volumeData = VolumeType::New();
  imgMsg = igtl::ImageMessage::New();

}

int Volumes::ITKtoIGTVolume()
{

  // Retrieve the image data from ITK image
  VolumeType::SizeType size;
  VolumeType::IndexType start;
  VolumeType::PointType origin;
  VolumeType::SpacingType spacing;

  size = this->volumeData->GetLargestPossibleRegion().GetSize();
  start = this->volumeData->GetLargestPossibleRegion().GetIndex();
  spacing = this->volumeData->GetSpacing();
  origin = this->volumeData->GetOrigin();

  // Set volume data to image message
  int   sizeIGT[3];
  int   scalarType;
  float originIGT[3];
  float spacingIGT[3];    // spacing (mm/pixel)
  //int   svsize[3];      // sub-volume size
  //int   svoffset[3];    // sub-volume offset
  sizeIGT[0] = size[0];         sizeIGT[1] = size[1];       sizeIGT[2] = size[2];
  spacingIGT[0] = spacing[0];   spacingIGT[1] = spacing[1]; spacingIGT[2] = spacing[2];
  originIGT[0] = origin[0]+((size[0]-1)*spacing[0]/2);      originIGT[1] = origin[1]+((size[1]-1)*spacing[1]/2);    originIGT[2] = origin[2]+((size[2]-1)*spacing[2]/2);
  //svsize[0] = sizeIGT[0];     svsize[1] = sizeIGT[1];     svsize[2] = sizeIGT[2];
  //svoffset[0] = start[0];     svoffset[1] = start[1];     svoffset[2] = start[2];
  scalarType = igtl::ImageMessage::TYPE_UINT8;

  this->imgMsg->SetDimensions( sizeIGT );
  this->imgMsg->SetSpacing( spacingIGT );
  this->imgMsg->SetScalarType( scalarType );
  this->imgMsg->SetOrigin( originIGT );
  this->imgMsg->AllocateScalars();

  // Set copy image data into ITK image
  memcpy( this->imgMsg->GetScalarPointer(), volumeData->GetBufferPointer(), this->imgMsg->GetImageSize() );
  // Pack (serialize) and send
  this->imgMsg->Pack();

  return 1;

}

class Clients
{

  public:

  Clients( char*, int );

  igtl::ClientSocket::Pointer socket;

};

Clients::Clients( char* host, int port )
{

  //Open a socket for th client
  socket = igtl::ClientSocket::New();

  //Connect to the server
  int r = socket->ConnectToServer( host, port );

  // Check the connection
  if ( (r) != 0 )
  {
    std::cerr << "Cannot connect to the server." << std::endl;
    exit(0);
  }

  std::cerr << "Client is connected to server: " << host << ":" << port << std::endl;

}

int RegisteredImage( Images* movingImagep, Images* secondImagep, Images* registeredImagep, RegistrationType::Pointer registration )
{

  // Use resulting transform from the registration to map the moving image into the moving/fixed image space
  ResampleFilterType::Pointer   resampler = ResampleFilterType::New();

  // Set moving image as input
  resampler->SetInput( movingImagep->imageData );

  // The Transform produced by the Registration method is passed into the resampling filter
  resampler->SetTransform( registration->GetOutput()->Get() );

  // Specifying parameters of the output image (default pixel value is set to gray in order to highlight the regions that are mapped outside of the moving image)
  resampler->SetSize( secondImagep->imageData->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin( secondImagep->imageData->GetOrigin() );
  resampler->SetOutputSpacing( secondImagep->imageData->GetSpacing() );
  resampler->SetOutputDirection( secondImagep->imageData->GetDirection() );
  resampler->SetDefaultPixelValue( 0 );
  resampler->Update();

  // Create registered ITKimage
  registeredImagep->imageData = resampler->GetOutput();

  // Set orientation of the registered image
  igtl::Matrix4x4           matrix;
  secondImagep->imgMsg->GetMatrix( matrix );
  registeredImagep->imgMsg->SetMatrix( matrix );

  return 1;

}

double RegistrationFunction( Images* fixedImagep, Images* movingImagep, Images* registeredImagep, ParametersType* finalParameters, int last )
{

  // Create components registration function
  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  // Each component is now connected to the instance of the registration method
  registration->SetMetric( metric );
  registration->SetOptimizer( optimizer );
  registration->SetTransform( transform );
  registration->SetInterpolator( interpolator );

  // Set the registration inputs
  registration->SetFixedImage( fixedImagep->imageData );
  registration->SetMovingImage( movingImagep->imageData );
  registration->SetFixedImageRegion( fixedImagep->imageData->GetLargestPossibleRegion() );

  //  Initialize the transform
  ParametersType initialParameters( transform->GetNumberOfParameters() );

  // translation matrix
  initialParameters[0] = 0.0;  // R(0,0)
  initialParameters[1] = 0.0;  // R(0,1)
  // No rotation allowed in this registration
  //initialParameters[2] = 0.0;
  //initialParameters[3] = 1.0;
  //initialParameters[4] = 0.0;
  //initialParameters[5] = 0.0;

  registration->SetInitialTransformParameters( initialParameters );

  // Set step size criteria
  optimizer->SetMaximumStepLength( 0.2 ); // If this is set too high, you will get a "itk::ERROR: MeanSquaresImageToImageMetric(0xa27ce70): Too many samples map outside moving image buffer: 1818 / 10000" error
  optimizer->SetMinimumStepLength( 0.05 );

  // Set a stopping criterion
  optimizer->SetNumberOfIterations( 100 );

  // Do the registration
  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 0;
  }

  // Create registered version of moving image
  if ( last == 1 )
  {
    if ( 0 == RegisteredImage( movingImagep, fixedImagep, registeredImagep, registration ) )
    {
      std::cerr << "Registering Image failed" << std::endl;
    }
  }

  //  The result of the registration process is an array of parameters that defines the spatial transformation in an unique way
  *finalParameters = registration->GetLastTransformParameters();
  std::cout << "Final parameters 2D/2D-Registration: " << finalParameters << std::endl;

  //  The value of the image metric corresponding to the last set of parameters can be obtained with the \code{GetValue()} method of the optimizer.
  std::cout << "Metric value: " << optimizer->GetValue() << std::endl;

  // Return matric value
  double val = (double)optimizer->GetValue();
  return val;

}

int LoadVolume( char* filename, Volumes* volume )
{

  // Open file reader and set NRRD file as input
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(filename);

  try
  {
    // Read volume from file
    reader->Update();
    volume->volumeData = reader->GetOutput();
    std::cerr << "Volume dimensions: " << volume->volumeData->GetLargestPossibleRegion().GetSize() << std::endl;
  }
  catch (itk::ExceptionObject &ex)
  {
    std::cerr << ex << std::endl;
    return 0;
  }

  std::cerr << "Volume Loaded" << std::endl;

  return 1;

}

int resliceImageVolume( Volumes* volume, int dStart[3], int dsize[2], Images* sliceImage )
{

  // Set parameters for desired image (reslice/ crop)
  VolumeType::IndexType         desiredStart;
  VolumeType::SizeType          desiredSize ;
  desiredStart[0] = dStart[0];  desiredStart[1] = dStart[1];    desiredStart[2] = dStart[2];
  desiredSize[0] = dsize[0];    desiredSize[1] = dsize[1];      desiredSize[2] = 0;

  VolumeType::RegionType        desiredRegion( desiredStart, desiredSize );
  std::cout << "Desired Region: " << desiredRegion << std::endl;
  std::cout << "Desired Region2: " << dStart[1] << std::endl;
  std::cout << "Desired Region2: " << desiredStart[1] << std::endl;

  // Create cropping/reslice
  FilterType::Pointer           filter = FilterType::New();
  filter->SetExtractionRegion( desiredRegion );
  filter->SetInput( volume->volumeData );
  filter->SetDirectionCollapseToIdentity(); // This is required.
  filter->Update();

  // Set resliced image to ITK image
  sliceImage->imageData = filter->GetOutput();

  // Get and set parameters for reslices image
  VolumeType::SpacingType       spacing;       // spacing (mm/pixel)
  float                         spacingsliceImage[3];
  VolumeType::IndexType         start1;
  VolumeType::PointType         origin;
  float                         startsliceImage[3];
  float                         originsliceImage[3];
  VolumeType::SizeType          size1;

  start1 = volume->volumeData->GetLargestPossibleRegion().GetIndex();
  spacing = volume->volumeData->GetSpacing();
  origin = volume->volumeData->GetOrigin();
  size1 = volume->volumeData->GetLargestPossibleRegion().GetSize();
  startsliceImage[0] = start1[0]+dStart[0]; startsliceImage[1] = start1[1]+dStart[1];   startsliceImage[2] = start1[2]+dStart[2];
  spacingsliceImage[0] = spacing[0];        spacingsliceImage[1] = spacing[1];          spacingsliceImage[2] = spacing[2];
  originsliceImage[0] = origin[0] + ( (size1[0]-1) * spacing[0]/2 ) + dStart[0];         originsliceImage[1] = origin[1] + ( (size1[1]-1) * spacing[1]/2 ) + dStart[1];       originsliceImage[2] = origin[2] + ( (size1[2]-1) * spacing[2]/2 ) + dStart[2];

  sliceImage->imageData->SetSpacing( spacingsliceImage );
  sliceImage->imageData->SetOrigin( originsliceImage );
  sliceImage->imgMsg->SetOrigin( originsliceImage );
  sliceImage->imgMsg->SetSpacing( spacingsliceImage );

  return 1;

}

int main(int argc, char* argv[])
{
  // Measure running time
  const clock_t begin_time = clock();

  if (argc != 5) // check number of arguments
  {
    // If not correct, print usage
    std::cerr << "    <filenameVolume>    : Filename of Volume.nrrd"            << std::endl;
    std::cerr << "    <hostnameSender1>   : IP or host name"                    << std::endl;
    std::cerr << "    <portSender1>       : Portimage # (18944 default)"        << std::endl;
    std::cerr << "    <sliceNumber>       : Number between 5 and 20 for now"    << std::endl;
    exit(0);
  }

  // Load volume data (NRRD file)
  char*     file = argv[1];
  Volumes   volume;
  LoadVolume( file, &volume );

  // Establish connections with server
  Clients   client1( argv[2], atoi(argv[3]) );

  // Send volume to slicer
  volume.ITKtoIGTVolume();
  client1.socket->Send( volume.imgMsg->GetPackPointer(), volume.imgMsg->GetPackSize() );

  // Create images for registration
  Images    fixedImage;
  Images    movingImage;
  Images    registeredImage;

  // Set parameters for testing resliceImage volume, get (part of) a slice from the volume
  int       sliceNumber = atoi(argv[4]);
  int       translation[2];     translation[0] = 5;                 translation[1] = 0; // Optional
  int       dStart[3];          dStart[0] = translation[0];         dStart[1] = translation[1];        dStart[2] =sliceNumber; // slice number "slicenumber" without the 5 most left pixels (so translated to the left)
  VolumeType::SizeType          size = volume.volumeData->GetLargestPossibleRegion().GetSize();
  int       dSize[2];           dSize[0] = size[0]-translation[0];  dSize[1] = size[1]-translation[1]; // Note the switch in axis for the translation

  // Get (part of) a slice of the volume (by cropping) and set it as fixed image for registration
  resliceImageVolume( &volume, dStart, dSize, &fixedImage );

  // Send the fixed image to Slicer
  fixedImage.ITKtoIGTImage();
  fixedImage.imgMsg->SetDeviceName( "imageSlice" );
  fixedImage.imgMsg->Pack();
  client1.socket->Send( fixedImage.imgMsg->GetPackPointer(), fixedImage.imgMsg->GetPackSize() );

  // Start registration of the resliced images with other resliced images in the area around it (5 slices before and 5 after) to find the best match
  int               testNumber = 11;
  double            values[testNumber];
  double            bestMatch[2];
  ParametersType    finalParameters;
  for ( int i = 0; i < testNumber; i++ )
  {
    // Measure running time for one reslice and registration loop
    const clock_t begin_time_reg = clock();

    // Set parameters for the reslice to match with the fixed image
    dStart[0] = 0;  dStart[1] = 0;  dStart[2] = sliceNumber - ((testNumber-1)/2)+i;

    //VolumeType::SizeType size = volume.volumeData->GetLargestPossibleRegion().GetSize();
    //dSize[2]; dSize[0]=size[0]; dSize[1]=size[1];

    // Get the new reslice of the volume (by cropping) and set it as the moving image for registration
    resliceImageVolume( &volume, dStart, dSize, &movingImage );

    // Register the moving and the fixed image and save th emetric values to find the best match
    values[i] = RegistrationFunction( &fixedImage, &movingImage, &registeredImage, &finalParameters, 0 );
    if ( i > 0 )
    {
      if( values[i] < bestMatch[1] )
      {
        bestMatch[0] = i;
        bestMatch[1] = values[i];
        std::cerr << "Best value: " << bestMatch[1] << std::endl;
      }
    }
    else
    {
      bestMatch[0] = i;
      bestMatch[1] = values[i];
      std::cerr << "Best value: " << bestMatch[1] << std::endl;
    }
    std::cerr << "Registration " << i+1 << " is done." << std::endl;
    std::cout << "Registration time: " << float( clock () - begin_time_reg) / CLOCKS_PER_SEC << std::endl;
  }

  // Give the best mathcing slice and the associated metric value and reistration values
  dStart[0] = 0;    dStart[1] = 0;    dStart[2] = sliceNumber - ((testNumber-1)/2) + bestMatch[0];
  resliceImageVolume( &volume, dStart, dSize, &movingImage );
  RegistrationFunction( &fixedImage, &movingImage, &registeredImage, &finalParameters, 1 ); // the 1 activated the image registration function
  std::cerr<< "Slice number " << sliceNumber << " matches best with slice number " << sliceNumber-((testNumber-1)/2)+bestMatch[0] << " with metric: " << bestMatch[1] << std::endl;
  std::cerr<< "Set translation: " << translation[0] << ", " << translation[1] << " vs. found translation by registration: " << finalParameters[0] << ", " << finalParameters[1] << std::endl;

  // Send the unregistered best matching slice and the registered best slice to Slicer
  movingImage.ITKtoIGTImage();
  movingImage.imgMsg->SetDeviceName( "bestMatch" );
  movingImage.imgMsg->Pack();
  client1.socket->Send( movingImage.imgMsg->GetPackPointer(), movingImage.imgMsg->GetPackSize() );
  registeredImage.ITKtoIGTImage();
  registeredImage.imgMsg->SetDeviceName( "registeredImage" );
  registeredImage.imgMsg->Pack();
  client1.socket->Send( registeredImage.imgMsg->GetPackPointer(), registeredImage.imgMsg->GetPackSize() );

  // Close connection
  client1.socket->CloseSocket();

  std::cout << "Total running time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;

}

