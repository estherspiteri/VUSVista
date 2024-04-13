import { IFile, ISample, ISampleVariant } from "../../models/view-samples.model";
import {
  IAddAcmgRuleRequest,
  IAddAcmgRuleResponse,
  IAddPhenotypeRequest,
  IAddPhenotypeResponse,
  IGetHPOTermsRequest,
  IGetHPOTermsResponse,
  IGetSampleRequest,
  IGetSampleResponse,
  ILoadAllSamplesResponse,
  IRemoveAcmgRuleRequest,
  IRemoveAcmgRuleResponse,
  IRemovePhenotypeRequest,
  IRemovePhenotypeResponse,
} from "./sample.dto";

export class SampleService {
  async loadAllSamples(): Promise<ILoadAllSamplesResponse> {
    const result: ILoadAllSamplesResponse = await fetch(`/sample/view`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getSample(input: IGetSampleRequest): Promise<IGetSampleResponse> {
    const result: IGetSampleResponse = await fetch(
      `/sample/view/${input.sampleId}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    let variants: ISampleVariant[] = []

    result.sample.variants.forEach((v) => {
      let files: IFile[] = [];

      v.files.forEach((f) => {
        files = files.concat({
          ...f,
          dateOfFileUpload: new Date(f.dateOfFileUpload),
        });
      });

      variants = variants.concat({...v, files: files})
    });

    let updatedSample: ISample = {
      ...result.sample,
      variants: variants,
    };

    return {
      isSuccess: result.isSuccess,
      sample: updatedSample,
      acmgRules: result.acmgRules,
    };
  }

  async addAcmgRule(input: IAddAcmgRuleRequest): Promise<IAddAcmgRuleResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleId", input.sampleId);
    data.append("variantId", input.variantId.toString());
    data.append("ruleId", input.ruleId.toString());

    const result: IAddAcmgRuleResponse = await fetch(`/sample/add-acmg-rule`, {
      method: "POST",
      body: data,
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async removeAcmgRule(
    input: IRemoveAcmgRuleRequest
  ): Promise<IRemoveAcmgRuleResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleId", input.sampleId);
    data.append("variantId", input.variantId.toString());
    data.append("ruleId", input.ruleId.toString());

    const result: IRemoveAcmgRuleResponse = await fetch(
      `/sample/remove-acmg-rule`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async getHPOTerms(input: IGetHPOTermsRequest): Promise<IGetHPOTermsResponse> {
    const result: IGetHPOTermsResponse = await fetch(
      `/sample/phenotype/${input.phenotype}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async addPhenotype(
    input: IAddPhenotypeRequest
  ): Promise<IAddPhenotypeResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleId", input.sampleId);
    data.append("phenotype", JSON.stringify(input.phenotype));

    const result: IGetHPOTermsResponse = await fetch(`/sample/add-phenotype`, {
      method: "POST",
      body: data,
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async removePhenotype(
    input: IRemovePhenotypeRequest
  ): Promise<IRemovePhenotypeResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleId", input.sampleId);
    data.append("phenotype", JSON.stringify(input.phenotype));

    const result: IGetHPOTermsResponse = await fetch(
      `/sample/remove-phenotype`,
      {
        method: "POST",
        body: data,
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const samplesService = new SampleService();
