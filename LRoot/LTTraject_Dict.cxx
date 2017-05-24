// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LTTraject_Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "LTTraject.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_LTTraject(void *p = 0);
   static void *newArray_LTTraject(Long_t size, void *p);
   static void delete_LTTraject(void *p);
   static void deleteArray_LTTraject(void *p);
   static void destruct_LTTraject(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LTTraject*)
   {
      ::LTTraject *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LTTraject >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LTTraject", ::LTTraject::Class_Version(), "LTTraject.hxx", 20,
                  typeid(::LTTraject), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LTTraject::Dictionary, isa_proxy, 4,
                  sizeof(::LTTraject) );
      instance.SetNew(&new_LTTraject);
      instance.SetNewArray(&newArray_LTTraject);
      instance.SetDelete(&delete_LTTraject);
      instance.SetDeleteArray(&deleteArray_LTTraject);
      instance.SetDestructor(&destruct_LTTraject);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LTTraject*)
   {
      return GenerateInitInstanceLocal((::LTTraject*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::LTTraject*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr LTTraject::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LTTraject::Class_Name()
{
   return "LTTraject";
}

//______________________________________________________________________________
const char *LTTraject::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTTraject*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LTTraject::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTTraject*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LTTraject::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTTraject*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LTTraject::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTTraject*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void LTTraject::Streamer(TBuffer &R__b)
{
   // Stream an object of class LTTraject.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LTTraject::Class(),this);
   } else {
      R__b.WriteClassBuffer(LTTraject::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LTTraject(void *p) {
      return  p ? new(p) ::LTTraject : new ::LTTraject;
   }
   static void *newArray_LTTraject(Long_t nElements, void *p) {
      return p ? new(p) ::LTTraject[nElements] : new ::LTTraject[nElements];
   }
   // Wrapper around operator delete
   static void delete_LTTraject(void *p) {
      delete ((::LTTraject*)p);
   }
   static void deleteArray_LTTraject(void *p) {
      delete [] ((::LTTraject*)p);
   }
   static void destruct_LTTraject(void *p) {
      typedef ::LTTraject current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LTTraject

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 210,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlETVector3gR_Dictionary();
   static void vectorlETVector3gR_TClassManip(TClass*);
   static void *new_vectorlETVector3gR(void *p = 0);
   static void *newArray_vectorlETVector3gR(Long_t size, void *p);
   static void delete_vectorlETVector3gR(void *p);
   static void deleteArray_vectorlETVector3gR(void *p);
   static void destruct_vectorlETVector3gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TVector3>*)
   {
      vector<TVector3> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TVector3>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TVector3>", -2, "vector", 210,
                  typeid(vector<TVector3>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETVector3gR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TVector3>) );
      instance.SetNew(&new_vectorlETVector3gR);
      instance.SetNewArray(&newArray_vectorlETVector3gR);
      instance.SetDelete(&delete_vectorlETVector3gR);
      instance.SetDeleteArray(&deleteArray_vectorlETVector3gR);
      instance.SetDestructor(&destruct_vectorlETVector3gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TVector3> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<TVector3>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETVector3gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TVector3>*)0x0)->GetClass();
      vectorlETVector3gR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETVector3gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETVector3gR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TVector3> : new vector<TVector3>;
   }
   static void *newArray_vectorlETVector3gR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TVector3>[nElements] : new vector<TVector3>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETVector3gR(void *p) {
      delete ((vector<TVector3>*)p);
   }
   static void deleteArray_vectorlETVector3gR(void *p) {
      delete [] ((vector<TVector3>*)p);
   }
   static void destruct_vectorlETVector3gR(void *p) {
      typedef vector<TVector3> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TVector3>

namespace ROOT {
   static TClass *vectorlETLorentzVectorgR_Dictionary();
   static void vectorlETLorentzVectorgR_TClassManip(TClass*);
   static void *new_vectorlETLorentzVectorgR(void *p = 0);
   static void *newArray_vectorlETLorentzVectorgR(Long_t size, void *p);
   static void delete_vectorlETLorentzVectorgR(void *p);
   static void deleteArray_vectorlETLorentzVectorgR(void *p);
   static void destruct_vectorlETLorentzVectorgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TLorentzVector>*)
   {
      vector<TLorentzVector> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TLorentzVector>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TLorentzVector>", -2, "vector", 210,
                  typeid(vector<TLorentzVector>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETLorentzVectorgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TLorentzVector>) );
      instance.SetNew(&new_vectorlETLorentzVectorgR);
      instance.SetNewArray(&newArray_vectorlETLorentzVectorgR);
      instance.SetDelete(&delete_vectorlETLorentzVectorgR);
      instance.SetDeleteArray(&deleteArray_vectorlETLorentzVectorgR);
      instance.SetDestructor(&destruct_vectorlETLorentzVectorgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TLorentzVector> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<TLorentzVector>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETLorentzVectorgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TLorentzVector>*)0x0)->GetClass();
      vectorlETLorentzVectorgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETLorentzVectorgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETLorentzVectorgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TLorentzVector> : new vector<TLorentzVector>;
   }
   static void *newArray_vectorlETLorentzVectorgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TLorentzVector>[nElements] : new vector<TLorentzVector>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETLorentzVectorgR(void *p) {
      delete ((vector<TLorentzVector>*)p);
   }
   static void deleteArray_vectorlETLorentzVectorgR(void *p) {
      delete [] ((vector<TLorentzVector>*)p);
   }
   static void destruct_vectorlETLorentzVectorgR(void *p) {
      typedef vector<TLorentzVector> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TLorentzVector>

namespace ROOT {
   static TClass *vectorlELTHitSegmentgR_Dictionary();
   static void vectorlELTHitSegmentgR_TClassManip(TClass*);
   static void *new_vectorlELTHitSegmentgR(void *p = 0);
   static void *newArray_vectorlELTHitSegmentgR(Long_t size, void *p);
   static void delete_vectorlELTHitSegmentgR(void *p);
   static void deleteArray_vectorlELTHitSegmentgR(void *p);
   static void destruct_vectorlELTHitSegmentgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LTHitSegment>*)
   {
      vector<LTHitSegment> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LTHitSegment>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LTHitSegment>", -2, "vector", 210,
                  typeid(vector<LTHitSegment>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELTHitSegmentgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<LTHitSegment>) );
      instance.SetNew(&new_vectorlELTHitSegmentgR);
      instance.SetNewArray(&newArray_vectorlELTHitSegmentgR);
      instance.SetDelete(&delete_vectorlELTHitSegmentgR);
      instance.SetDeleteArray(&deleteArray_vectorlELTHitSegmentgR);
      instance.SetDestructor(&destruct_vectorlELTHitSegmentgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LTHitSegment> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<LTHitSegment>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELTHitSegmentgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LTHitSegment>*)0x0)->GetClass();
      vectorlELTHitSegmentgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELTHitSegmentgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELTHitSegmentgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LTHitSegment> : new vector<LTHitSegment>;
   }
   static void *newArray_vectorlELTHitSegmentgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LTHitSegment>[nElements] : new vector<LTHitSegment>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELTHitSegmentgR(void *p) {
      delete ((vector<LTHitSegment>*)p);
   }
   static void deleteArray_vectorlELTHitSegmentgR(void *p) {
      delete [] ((vector<LTHitSegment>*)p);
   }
   static void destruct_vectorlELTHitSegmentgR(void *p) {
      typedef vector<LTHitSegment> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LTHitSegment>

namespace {
  void TriggerDictionaryInitialization_LTTraject_Dict_Impl() {
    static const char* headers[] = {
"LTTraject.hxx",
0
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/.",
"/usr/local/root-6.08.00-build/include",
"/home/nmcfadde/doubleBeta/LRoot/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LTTraject_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$LTTraject.hxx")))  TLorentzVector;
class __attribute__((annotate("$clingAutoload$LTTraject.hxx")))  TVector3;
class __attribute__((annotate("$clingAutoload$LTTraject.hxx")))  LTHitSegment;
class __attribute__((annotate("$clingAutoload$LTTraject.hxx")))  LTTraject;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LTTraject_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "LTTraject.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"LTTraject", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LTTraject_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LTTraject_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LTTraject_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LTTraject_Dict() {
  TriggerDictionaryInitialization_LTTraject_Dict_Impl();
}
