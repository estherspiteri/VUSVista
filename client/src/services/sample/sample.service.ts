import {
  IAddAcmgRuleRequest,
  IGetSampleRequest,
  IGetSampleResponse,
  ILoadAllSamplesResponse,
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

    let updatedSampleList = result?.sampleList?.map((s) => {
      return { ...s, dateOfFileUpload: new Date(s.dateOfFileUpload) };
    });

    return { isSuccess: result.isSuccess, sampleList: updatedSampleList };
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

    let updatedSample = {
      ...result.sample,
      dateOfFileUpload: new Date(result.sample.dateOfFileUpload),
    };

    return {
      isSuccess: result.isSuccess,
      sample: updatedSample,
      acmgRuleNames: result.acmgRuleNames,
    };
  }

  async addAcmgRule(input: IAddAcmgRuleRequest) {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleId", input.sampleId);
    data.append("variantId", input.variantId.toString());
    data.append("ruleName", input.ruleName);

    await fetch(`/sample/add-acmg-rule`, {
      method: "POST",
      body: data,
    }).catch((error) => console.error("error============:", error)); //TODO: handle error
  }

  async removeAcmgRule(input: IAddAcmgRuleRequest) {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("sampleId", input.sampleId);
    data.append("variantId", input.variantId.toString());
    data.append("ruleName", input.ruleName);

    await fetch(`/sample/remove-acmg-rule`, {
      method: "POST",
      body: data,
    }).catch((error) => console.error("error============:", error)); //TODO: handle error
  }
}

export const samplesService = new SampleService();
